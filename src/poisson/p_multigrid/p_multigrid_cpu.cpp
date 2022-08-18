#include "p_multigrid.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

#include "dg_utils.h"

// Copy PETSc vec array to OP2 dat
void PMultigrid::copy_vec_to_dat(op_dat dat, const double *dat_d) {
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_WRITE),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 2, copy_args);

  int setSize = dat->set->size;
  const int *p = (int *)mesh->order->data;

  int vec_ind = 0;
  int block_start = 0;
  int block_count = 0;
  for(int i = 0; i < setSize; i++) {
    const int N = p[i];

    if(N == DG_ORDER) {
      if(block_count == 0) {
        block_start = i;
        block_count++;
        continue;
      } else {
        block_count++;
        continue;
      }
    } else {
      if(block_count != 0) {
        double *block_start_dat_c = (double *)dat->data + block_start * dat->dim;
        memcpy(block_start_dat_c, dat_d + vec_ind, block_count * DG_NP * sizeof(double));
        vec_ind += DG_NP * block_count;
      }
      block_count = 0;
    }

    double *v_c = (double *)dat->data + i * dat->dim;
    int Np, Nfp;
    DGUtils::basic_constants(N, &Np, &Nfp);

    memcpy(v_c, dat_d + vec_ind, Np * sizeof(double));
    vec_ind += Np;
  }

  if(block_count != 0) {
    double *block_start_dat_c = (double *)dat->data + block_start * dat->dim;
    memcpy(block_start_dat_c, dat_d + vec_ind, block_count * DG_NP * sizeof(double));
    vec_ind += DG_NP * block_count;
  }

  op_mpi_set_dirtybit(2, copy_args);
}

// Copy OP2 dat to PETSc vec array
void PMultigrid::copy_dat_to_vec(op_dat dat, double *dat_d) {
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 2, copy_args);

  int setSize = dat->set->size;
  const int *p = (int *)mesh->order->data;

  int vec_ind = 0;
  int block_start = 0;
  int block_count = 0;
  for(int i = 0; i < setSize; i++) {
    const int N = p[i];

    if(N == DG_ORDER) {
      if(block_count == 0) {
        block_start = i;
        block_count++;
        continue;
      } else {
        block_count++;
        continue;
      }
    } else {
      if(block_count != 0) {
        const double *block_start_dat_c = (double *)dat->data + block_start * dat->dim;
        memcpy(dat_d + vec_ind, block_start_dat_c, block_count * DG_NP * sizeof(double));
        vec_ind += DG_NP * block_count;
      }
      block_count = 0;
    }

    const double *v_c = (double *)dat->data + i * dat->dim;
    int Np, Nfp;
    DGUtils::basic_constants(N, &Np, &Nfp);

    memcpy(dat_d + vec_ind, v_c, Np * sizeof(double));
    vec_ind += Np;
  }

  if(block_count != 0) {
    const double *block_start_dat_c = (double *)dat->data + block_start * dat->dim;
    memcpy(dat_d + vec_ind, block_start_dat_c, block_count * DG_NP * sizeof(double));
    vec_ind += DG_NP * block_count;
  }

  op_mpi_set_dirtybit(2, copy_args);
}

// Create a PETSc vector for CPUs
void PMultigrid::create_vec(Vec *v) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, pMatrix->unknowns, PETSC_DECIDE);
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void PMultigrid::load_vec(Vec *v, op_dat v_dat) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);

  copy_dat_to_vec(v_dat, v_ptr);

  VecRestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void PMultigrid::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);

  copy_vec_to_dat(v_dat, v_ptr);

  VecRestoreArrayRead(*v, &v_ptr);
}


PetscErrorCode matBMult(Mat A, Vec x, Vec y) {
  PMultigrid *p_multigrid;
  MatShellGetContext(A, &p_multigrid);
  const double *x_ptr;
  double *y_ptr;
  VecGetArrayRead(x, &x_ptr);
  VecGetArray(y, &y_ptr);

  p_multigrid->calc_rhs(x_ptr, y_ptr);

  VecRestoreArrayRead(x, &x_ptr);
  VecRestoreArray(y, &y_ptr);
  return 0;
}

void PMultigrid::create_shell_mat(Mat *mat) {
  MatCreateShell(PETSC_COMM_WORLD, pMatrix->unknowns, pMatrix->unknowns, PETSC_DETERMINE, PETSC_DETERMINE, this, mat);
  MatShellSetOperation(*mat, MATOP_MULT, (void(*)(void))matBMult);
  MatShellSetVecType(*mat, VECSTANDARD);
}

void PMultigrid::setMatrix() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreate(PETSC_COMM_WORLD, &pMat);
  pMatInit = true;
  MatSetSizes(pMat, pMatrix->unknowns, pMatrix->unknowns, PETSC_DECIDE, PETSC_DECIDE);

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
    op_arg_dat(pMatrix->glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 3, args);
  const double *op1_data = (double *)pMatrix->op1->data;
  const int *glb = (int *)pMatrix->glb_ind->data;
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
    op_arg_dat(pMatrix->glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(pMatrix->glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(pMatrix->orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(pMatrix->orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->edges, 6, edge_args);

  const double *op2L_data = (double *)pMatrix->op2[0]->data;
  const double *op2R_data = (double *)pMatrix->op2[1]->data;
  const int *glb_l = (int *)pMatrix->glb_indL->data;
  const int *glb_r = (int *)pMatrix->glb_indR->data;
  const int *p_l = (int *)pMatrix->orderL->data;
  const int *p_r = (int *)pMatrix->orderR->data;

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