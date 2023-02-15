#include "matrices/poisson_matrix.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

#include "dg_utils.h"
#include "dg_global_constants/dg_global_constants_2d.h"
#include "timing.h"

extern Timing *timer;

void get_num_nodes(const int N, int *Np, int *Nfp) {
  #if DG_DIM == 2
  DGUtils::numNodes2D(N, Np, Nfp);
  #elif DG_DIM == 3
  DGUtils::numNodes3D(N, Np, Nfp);
  #endif
}

void PoissonMatrix::set_glb_ind() {
  timer->startTimer("PoissonMat - glb_ind");
  op_arg op2_args[] = {
    op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(_mesh->order->set, 1, op2_args);
  const int setSize_ = _mesh->order->set->size;
  int *tempOrder_ = (int *)malloc(setSize_ * sizeof(int));
  cudaMemcpy(tempOrder_, _mesh->order->data_d, setSize_ * sizeof(int), cudaMemcpyDeviceToHost);
  unknowns = 0;
  for(int i = 0; i < setSize_; i++) {
    int Np, Nfp;
    get_num_nodes(tempOrder_[i], &Np, &Nfp);
    unknowns += Np;
  }
  free(tempOrder_);
  op_mpi_set_dirtybit_cuda(1, op2_args);
  int global_ind = 0;
  #ifdef INS_MPI
  global_ind = get_global_mat_start_ind(unknowns);
  #endif
  op_arg args[] = {
    op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(_mesh->cells, 2, args);

  const int setSize = _mesh->cells->size;
  int *tempOrder = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(tempOrder, _mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  int *data_ptr = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(data_ptr, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  int ind = global_ind;
  for(int i = 0; i < _mesh->cells->size; i++) {
    int Np, Nfp;
    get_num_nodes(tempOrder[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  cudaMemcpy(glb_ind->data_d, data_ptr, setSize * sizeof(int), cudaMemcpyHostToDevice);

  op_mpi_set_dirtybit_cuda(2, args);
  free(data_ptr);
  free(tempOrder);
  timer->endTimer("PoissonMat - glb_ind");
}

void PoissonMatrix::setPETScMatrix() {
  if(!petscMatInit) {
    MatCreate(PETSC_COMM_WORLD, &pMat);
    petscMatInit = true;
    MatSetSizes(pMat, unknowns, unknowns, PETSC_DECIDE, PETSC_DECIDE);

    #ifdef INS_MPI
    MatSetType(pMat, MATMPIAIJCUSPARSE);
    MatMPIAIJSetPreallocation(pMat, DG_NP * (DG_NUM_FACES + 1), NULL, 0, NULL);
    #else
    MatSetType(pMat, MATSEQAIJCUSPARSE);
    MatSeqAIJSetPreallocation(pMat, DG_NP * (DG_NUM_FACES + 1), NULL);
    #endif
    MatSetOption(pMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  }
  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(_mesh->cells, 3, args);

  const int setSize = _mesh->cells->size;
  DG_FP *op1_data = (DG_FP *)malloc(DG_NP * DG_NP * setSize * sizeof(DG_FP));
  int *glb   = (int *)malloc(setSize * sizeof(int));
  int *order = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(op1_data, op1->data_d, setSize * DG_NP * DG_NP * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order, _mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  op_mpi_set_dirtybit_cuda(3, args);

  #ifdef DG_COL_MAJ
  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_FALSE);
  #else
  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_TRUE);
  #endif

  for(int i = 0; i < setSize; i++) {
    int Np, Nfp;
    get_num_nodes(order[i], &Np, &Nfp);
    int currentRow = glb[i];
    int currentCol = glb[i];
    int idxm[DG_NP], idxn[DG_NP];
    for(int n = 0; n < DG_NP; n++) {
      idxm[n] = currentRow + n;
      idxn[n] = currentCol + n;
    }

    MatSetValues(pMat, Np, idxm, Np, idxn, &op1_data[i * DG_NP * DG_NP], INSERT_VALUES);
  }

  free(op1_data);
  free(glb);
  free(order);

  op_arg edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(_mesh->faces, 6, edge_args);
  DG_FP *op2L_data = (DG_FP *)malloc(DG_NP * DG_NP * _mesh->faces->size * sizeof(DG_FP));
  DG_FP *op2R_data = (DG_FP *)malloc(DG_NP * DG_NP * _mesh->faces->size * sizeof(DG_FP));
  int *glb_l = (int *)malloc(_mesh->faces->size * sizeof(int));
  int *glb_r = (int *)malloc(_mesh->faces->size * sizeof(int));
  int *order_l = (int *)malloc(_mesh->faces->size * sizeof(int));
  int *order_r = (int *)malloc(_mesh->faces->size * sizeof(int));

  cudaMemcpy(op2L_data, op2[0]->data_d, DG_NP * DG_NP * _mesh->faces->size * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaMemcpy(op2R_data, op2[1]->data_d, DG_NP * DG_NP * _mesh->faces->size * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_l, glb_indL->data_d, _mesh->faces->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, _mesh->faces->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_l, orderL->data_d, _mesh->faces->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_r, orderR->data_d, _mesh->faces->size * sizeof(int), cudaMemcpyDeviceToHost);

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < _mesh->faces->size; i++) {
    int leftRow = glb_l[i];
    int rightRow = glb_r[i];
    int NpL, NpR, Nfp;
    get_num_nodes(order_l[i], &NpL, &Nfp);
    get_num_nodes(order_r[i], &NpR, &Nfp);

    int idxl[DG_NP], idxr[DG_NP];
    for(int n = 0; n < DG_NP; n++) {
      idxl[n] = leftRow + n;
      idxr[n] = rightRow + n;
    }

    MatSetValues(pMat, NpL, idxl, NpR, idxr, &op2L_data[i * DG_NP * DG_NP], INSERT_VALUES);
    MatSetValues(pMat, NpR, idxr, NpL, idxl, &op2R_data[i * DG_NP * DG_NP], INSERT_VALUES);
  }

  free(op2L_data);
  free(op2R_data);
  free(glb_l);
  free(glb_r);
  free(order_l);
  free(order_r);

  op_mpi_set_dirtybit_cuda(6, edge_args);

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
}
