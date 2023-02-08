#include "matrices/2d/poisson_matrix_2d.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

#include "dg_utils.h"
#include "dg_global_constants/dg_global_constants_2d.h"
#include "timing.h"

extern Timing *timer;

void PoissonMatrix2D::set_glb_ind() {
  int global_ind = 0;
  #ifdef INS_MPI
  timer->startTimer("PoissonMat - glb_ind");
  op_arg op2_args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->order->set, 1, op2_args);
  const int setSize = mesh->order->set->size;
  const int *tempOrder = (int *)mesh->order->data;
  unknowns = 0;
  for(int i = 0; i < setSize; i++) {
    int Np, Nfp;
    DGUtils::numNodes2D(tempOrder[i], &Np, &Nfp);
    unknowns += Np;
  }
  op_mpi_set_dirtybit(1, op2_args);
  global_ind = get_global_mat_start_ind(unknowns);
  #endif
  op_arg args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 2, args);

  const int *p = (int *)mesh->order->data;
  int *data_ptr = (int *)glb_ind->data;
  int ind = global_ind;
  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::numNodes2D(p[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  op_mpi_set_dirtybit(2, args);
  timer->endTimer("PoissonMat - glb_ind");
}

void PoissonMatrix2D::setPETScMatrix() {
  // TODO update for p-adpativity/p-multigrid
  const int dg_np = DG_CONSTANTS[(DG_ORDER - 1) * 5];
  if(!petscMatInit) {
    MatCreate(PETSC_COMM_WORLD, &pMat);
    petscMatInit = true;
    MatSetSizes(pMat, mesh->cells->size * dg_np, mesh->cells->size * dg_np, PETSC_DECIDE, PETSC_DECIDE);

    #ifdef INS_MPI
    MatSetType(pMat, MATMPIAIJ);
    MatMPIAIJSetPreallocation(pMat, dg_np * 4, NULL, 0, NULL);
    #else
    MatSetType(pMat, MATSEQAIJ);
    MatSeqAIJSetPreallocation(pMat, dg_np * 4, NULL);
    #endif
    MatSetOption(pMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  }

  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 3, args);
  const double *op1_data = (double *)op1->data;
  const int *glb = (int *)glb_ind->data;
  const int *p = (int *)mesh->order->data;

  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_FALSE);

  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::numNodes2D(p[i], &Np, &Nfp);
    int currentRow = glb[i];
    int currentCol = glb[i];

    int idxm[DG_NP], idxn[DG_NP];
    for(int n = 0; n < DG_NP; n++) {
      idxm[n] = currentRow + n;
      idxn[n] = currentCol + n;
    }

    MatSetValues(pMat, Np, idxm, Np, idxn, &op1_data[i * DG_NP * DG_NP], INSERT_VALUES);
  }

  op_mpi_set_dirtybit(3, args);

  op_arg edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->faces, 6, edge_args);

  const double *op2L_data = (double *)op2[0]->data;
  const double *op2R_data = (double *)op2[1]->data;
  const int *glb_l = (int *)glb_indL->data;
  const int *glb_r = (int *)glb_indR->data;
  const int *p_l = (int *)orderL->data;
  const int *p_r = (int *)orderR->data;

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < mesh->faces->size; i++) {
    int leftRow = glb_l[i];
    int rightRow = glb_r[i];
    int NpL, NpR, Nfp;
    DGUtils::numNodes2D(p_l[i], &NpL, &Nfp);
    DGUtils::numNodes2D(p_r[i], &NpR, &Nfp);

    int idxl[DG_NP], idxr[DG_NP];
    for(int n = 0; n < DG_NP; n++) {
      idxl[n] = leftRow + n;
      idxr[n] = rightRow + n;
    }

    MatSetValues(pMat, NpL, idxl, NpR, idxr, &op2L_data[i * DG_NP * DG_NP], INSERT_VALUES);
    MatSetValues(pMat, NpR, idxr, NpL, idxl, &op2R_data[i * DG_NP * DG_NP], INSERT_VALUES);
  }

  op_mpi_set_dirtybit(6, edge_args);

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
}
