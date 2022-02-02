#include "poisson.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

#include "dg_utils.h"

int PoissonSolve::get_local_unknowns() {
  op_arg op2_args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->order->set, 1, op2_args);
  const int setSize = mesh->order->set->size;
  const int *p = (int *)mesh->order->data;
  int local_unkowns = 0;
  for(int i = 0; i < setSize; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(p[i], &Np, &Nfp);
    local_unkowns += Np;
  }
  op_mpi_set_dirtybit(1, op2_args);
  return local_unkowns;
}

// Copy PETSc vec array to OP2 dat
void PoissonSolve::copy_vec_to_dat(op_dat dat, const double *dat_d) {
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
void PoissonSolve::copy_dat_to_vec(op_dat dat, double *dat_d) {
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
void PoissonSolve::create_vec(Vec *v) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, unknowns, PETSC_DECIDE);
}

// Destroy a PETSc vector
void PoissonSolve::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void PoissonSolve::load_vec(Vec *v, op_dat v_dat) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);

  copy_dat_to_vec(v_dat, v_ptr);

  VecRestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void PoissonSolve::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);

  copy_vec_to_dat(v_dat, v_ptr);

  VecRestoreArrayRead(*v, &v_ptr);
}

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
    op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->edges, 6, edge_args);

  const double *op2L_data = (double *)op2[0]->data;
  const double *op2R_data = (double *)op2[1]->data;
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

void PressureSolve::setAMGXMat() {

}

void PressureSolve::uploadAMGXVec(AMGX_vector_handle *vec, op_dat dat) {

}
