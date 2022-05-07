#include "poisson.h"
#include "poisson_HYPRE.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

#include "dg_utils.h"

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  PoissonSolve *poisson;
  MatShellGetContext(A, &poisson);
  const double *x_ptr;
  double *y_ptr;
  VecGetArrayRead(x, &x_ptr);
  VecGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecRestoreArrayRead(x, &x_ptr);
  VecRestoreArray(y, &y_ptr);
  return 0;
}

void PoissonSolve::create_shell_mat() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreateShell(PETSC_COMM_WORLD, unknowns, unknowns, PETSC_DETERMINE, PETSC_DETERMINE, this, &pMat);
  MatShellSetOperation(pMat, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(pMat, VECSTANDARD);

  pMatInit = true;
}

PetscErrorCode precon(PC pc, Vec x, Vec y) {
  PoissonSolve *poisson;
  PCShellGetContext(pc, (void **)&poisson);
  const double *x_ptr;
  double *y_ptr;
  VecGetArrayRead(x, &x_ptr);
  VecGetArray(y, &y_ptr);

  poisson->precond(x_ptr, y_ptr);

  VecRestoreArrayRead(x, &x_ptr);
  VecRestoreArray(y, &y_ptr);
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
  op_mpi_halo_exchanges(mesh->cells, 2, args);

  const int *p = (int *)mesh->order->data;
  int *data_ptr = (int *)glb_ind->data;
  int ind = global_ind;
  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(p[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  op_mpi_set_dirtybit(2, args);
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
  op_mpi_halo_exchanges(mesh->cells, 2, args);

  const int *p = (int *)mesh->order->data;
  int *data_ptr = (int *)glb_ind->data;
  int ind = global_ind;
  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(p[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  op_mpi_set_dirtybit(2, args);
}

void PoissonSolve::setMatrix() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreate(PETSC_COMM_WORLD, &pMat);
  pMatInit = true;
  MatSetSizes(pMat, unknowns, unknowns, PETSC_DECIDE, PETSC_DECIDE);

  #ifdef INS_MPI
  MatSetType(pMat, MATMPIAIJ);
  MatMPIAIJSetPreallocation(pMat, DG_NP * 4, NULL, 0, NULL);
  #else
  MatSetType(pMat, MATSEQAIJ);
  MatSeqAIJSetPreallocation(pMat, DG_NP * 4, NULL);
  #endif
  MatSetOption(pMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(pMatrix->op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 3, args);
  const double *op1_data = (double *)pMatrix->op1->data;
  const int *glb = (int *)glb_ind->data;
  const int *p = (int *)mesh->order->data;

  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_FALSE);

  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(p[i], &Np, &Nfp);
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
    op_arg_dat(pMatrix->op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(pMatrix->op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->edges, 6, edge_args);

  const double *op2L_data = (double *)pMatrix->op2[0]->data;
  const double *op2R_data = (double *)pMatrix->op2[1]->data;
  const int *glb_l = (int *)glb_indL->data;
  const int *glb_r = (int *)glb_indR->data;
  const int *p_l = (int *)orderL->data;
  const int *p_r = (int *)orderR->data;

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < mesh->edges->size; i++) {
    int leftRow = glb_l[i];
    int rightRow = glb_r[i];
    int NpL, NpR, Nfp;
    DGUtils::basic_constants(p_l[i], &NpL, &Nfp);
    DGUtils::basic_constants(p_r[i], &NpR, &Nfp);

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

void PoissonSolveHYPRE::setMatrix() {
  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(pMatrix->op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 3, args);

  const int setSize = mesh->cells->size;
  const int *glb    = (int *)glb_ind->data;
  const int *order  = (int *)mesh->order->data;

  int *idxm  = (int *)malloc(DG_NP * setSize * sizeof(int));
  int *idxn  = (int *)malloc(DG_NP * DG_NP * setSize * sizeof(int));
  int *nCols = (int *)malloc(DG_NP * setSize * sizeof(int));

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

  HYPRE_IJMatrixAddToValues(mat, DG_NP * setSize, nCols, idxm, idxn, (double *)pMatrix->op1->data);

  op_mpi_set_dirtybit(3, args);

  free(idxm);
  free(idxn);
  free(nCols);

  op_arg edge_args[] = {
    op_arg_dat(pMatrix->op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(pMatrix->op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->edges, 6, edge_args);
  const int *glb_l   = (int *)glb_indL->data;
  const int *glb_r   = (int *)glb_indR->data;
  const int *order_l = (int *)orderL->data;
  const int *order_r = (int *)orderR->data;

  int *idxlm  = (int *)malloc(DG_NP * mesh->edges->size * sizeof(int));
  int *idxln  = (int *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(int));
  int *idxrm  = (int *)malloc(DG_NP * mesh->edges->size * sizeof(int));
  int *idxrn  = (int *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(int));
  int *nCols2 = (int *)malloc(DG_NP * mesh->edges->size * sizeof(int));

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

  HYPRE_IJMatrixAddToValues(mat, DG_NP * mesh->edges->size, nCols2, idxlm, idxln, (double *)pMatrix->op2[0]->data);
  HYPRE_IJMatrixAddToValues(mat, DG_NP * mesh->edges->size, nCols2, idxrm, idxrn, (double *)pMatrix->op2[1]->data);

  free(idxlm);
  free(idxln);
  free(idxrm);
  free(idxrn);
  free(nCols2);

  op_mpi_set_dirtybit(6, edge_args);
}
