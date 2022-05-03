#include "poisson.h"
#include "poisson_HYPRE.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

#include "dg_utils.h"

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  PoissonSolve *poisson;
  MatShellGetContext(A, &poisson);
  const double *x_ptr;
  double *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);
  return 0;
}

void PoissonSolve::create_shell_mat() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreateShell(PETSC_COMM_WORLD, unknowns, unknowns, PETSC_DETERMINE, PETSC_DETERMINE, this, &pMat);
  MatShellSetOperation(pMat, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(pMat, VECCUDA);

  pMatInit = true;
}

PetscErrorCode precon(PC pc, Vec x, Vec y) {
  PoissonSolve *poisson;
  PCShellGetContext(pc, (void **)&poisson);
  const double *x_ptr;
  double *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  poisson->precond(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);
  return 0;
}

void PoissonSolve::set_shell_pc(PC pc) {
  PCShellSetApply(pc, precon);
  PCShellSetContext(pc, this);
}

void PoissonSolve::setGlbInd() {
  int global_ind = 0;
  #ifdef INS_MPI
  global_ind = get_global_mat_start_ind(unknowns);
  #endif
  op_arg args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 2, args);

  const int setSize = mesh->cells->size;
  int *tempOrder = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(tempOrder, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  int *data_ptr = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(data_ptr, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  int ind = global_ind;
  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(tempOrder[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  cudaMemcpy(glb_ind->data_d, data_ptr, setSize * sizeof(int), cudaMemcpyHostToDevice);

  op_mpi_set_dirtybit_cuda(2, args);
  free(data_ptr);
  free(tempOrder);
}

void PoissonSolveHYPRE::setGlbInd() {
  int global_ind = 0;
  #ifdef INS_MPI
  global_ind = get_global_mat_start_ind(unknowns);
  #endif
  op_arg args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 2, args);

  const int setSize = mesh->cells->size;
  int *tempOrder = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(tempOrder, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  int *data_ptr = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(data_ptr, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  int ind = global_ind;
  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(tempOrder[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  cudaMemcpy(glb_ind->data_d, data_ptr, setSize * sizeof(int), cudaMemcpyHostToDevice);

  op_mpi_set_dirtybit_cuda(2, args);
  free(data_ptr);
  free(tempOrder);
}

void PoissonSolve::setMatrix() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreate(PETSC_COMM_WORLD, &pMat);
  pMatInit = true;
  MatSetSizes(pMat, unknowns, unknowns, PETSC_DECIDE, PETSC_DECIDE);

  #ifdef INS_MPI
  MatSetType(pMat, MATMPIAIJCUSPARSE);
  MatMPIAIJSetPreallocation(pMat, DG_NP * 4, NULL, 0, NULL);
  #else
  MatSetType(pMat, MATSEQAIJCUSPARSE);
  MatSeqAIJSetPreallocation(pMat, DG_NP * 4, NULL);
  #endif
  MatSetOption(pMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(pMatrix->op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 3, args);

  const int setSize = mesh->cells->size;
  double *op1_data = (double *)malloc(DG_NP * DG_NP * setSize * sizeof(double));
  int *glb   = (int *)malloc(setSize * sizeof(int));
  int *order = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(op1_data, pMatrix->op1->data_d, setSize * DG_NP * DG_NP * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  op_mpi_set_dirtybit_cuda(3, args);

  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_FALSE);

  for(int i = 0; i < setSize; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(order[i], &Np, &Nfp);
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
    op_arg_dat(pMatrix->op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(pMatrix->op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->edges, 6, edge_args);
  double *op2L_data = (double *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(double));
  double *op2R_data = (double *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(double));
  int *glb_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *glb_r = (int *)malloc(mesh->edges->size * sizeof(int));
  int *order_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *order_r = (int *)malloc(mesh->edges->size * sizeof(int));

  cudaMemcpy(op2L_data, pMatrix->op2[0]->data_d, DG_NP * DG_NP * mesh->edges->size * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(op2R_data, pMatrix->op2[1]->data_d, DG_NP * DG_NP * mesh->edges->size * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_l, glb_indL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_l, orderL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_r, orderR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < mesh->edges->size; i++) {
    int leftRow = glb_l[i];
    int rightRow = glb_r[i];
    int NpL, NpR, Nfp;
    DGUtils::basic_constants(order_l[i], &NpL, &Nfp);
    DGUtils::basic_constants(order_r[i], &NpR, &Nfp);

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

void PoissonSolveHYPRE::setMatrix() {
  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(pMatrix->op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 3, args);

  const int setSize = mesh->cells->size;
  int *glb   = (int *)malloc(setSize * sizeof(int));
  int *order = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(glb, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  int *idxm, *idxn, *nCols;
  cudaMallocManaged(&idxm, DG_NP * setSize * sizeof(int));
  cudaMallocManaged(&idxn, DG_NP * DG_NP * setSize * sizeof(int));
  cudaMallocManaged(&nCols, DG_NP * setSize * sizeof(int));

  for(int i = 0; i < setSize; i++) {
    int currentRow = glb[i];
    int currentCol = glb[i];
    for(int n = 0; n < DG_NP; n++) {
      idxm[i * DG_NP + n] = currentRow + n;
      nCols[i * DG_NP + n] = DG_NP;
      for(int m = 0; m < DG_NP; m++) {
        idxn[i * DG_NP * DG_NP + n * DG_NP + m] = currentCol + m;
      }
    }
  }

  HYPRE_IJMatrixAddToValues(mat, DG_NP * setSize, nCols, idxm, idxn, (double *)pMatrix->op1->data_d);

  op_mpi_set_dirtybit_cuda(3, args);

  cudaFree(idxm);
  cudaFree(idxn);
  cudaFree(nCols);

  free(glb);
  free(order);

  op_arg edge_args[] = {
    op_arg_dat(pMatrix->op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(pMatrix->op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->edges, 6, edge_args);
  int *glb_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *glb_r = (int *)malloc(mesh->edges->size * sizeof(int));
  int *order_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *order_r = (int *)malloc(mesh->edges->size * sizeof(int));

  cudaMemcpy(glb_l, glb_indL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_l, orderL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_r, orderR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);

  int *idxlm, *idxln, *idxrm, *idxrn, *nCols2;
  cudaMallocManaged(&idxlm, DG_NP * mesh->edges->size * sizeof(int));
  cudaMallocManaged(&idxln, DG_NP * DG_NP * mesh->edges->size * sizeof(int));
  cudaMallocManaged(&idxrm, DG_NP * mesh->edges->size * sizeof(int));
  cudaMallocManaged(&idxrn, DG_NP * DG_NP * mesh->edges->size * sizeof(int));
  cudaMallocManaged(&nCols2, DG_NP * mesh->edges->size * sizeof(int));

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < mesh->edges->size; i++) {
    int leftRow = glb_l[i];
    int rightRow = glb_r[i];

    for(int n = 0; n < DG_NP; n++) {
      idxlm[i * DG_NP + n] = leftRow + n;
      idxrm[i * DG_NP + n] = rightRow + n;
      nCols2[i * DG_NP + n] = DG_NP;
      for(int m = 0; m < DG_NP; m++) {
        idxln[i * DG_NP * DG_NP + n * DG_NP + m] = rightRow + m;
        idxrn[i * DG_NP * DG_NP + n * DG_NP + m] = leftRow + m;
      }
    }
  }

  HYPRE_IJMatrixAddToValues(mat, DG_NP * mesh->edges->size, nCols2, idxlm, idxln, (double *)pMatrix->op2[0]->data_d);
  HYPRE_IJMatrixAddToValues(mat, DG_NP * mesh->edges->size, nCols2, idxrm, idxrn, (double *)pMatrix->op2[1]->data_d);

  cudaFree(idxlm);
  cudaFree(idxln);
  cudaFree(idxrm);
  cudaFree(idxrn);
  cudaFree(nCols2);

  free(glb_l);
  free(glb_r);
  free(order_l);
  free(order_r);

  op_mpi_set_dirtybit_cuda(6, edge_args);
}
