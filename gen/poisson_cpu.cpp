#include "poisson.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

#include "dg_compiler_defs.h"

// Copy PETSc vec array to OP2 dat
void PoissonSolve::copy_vec_to_dat(op_dat dat, const double *dat_d) {
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, 10, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(dat->set, 1, copy_args);
  memcpy(dat->data, dat_d, dat->set->size * 10 * sizeof(double));
  op_mpi_set_dirtybit(1, copy_args);
}

// Copy OP2 dat to PETSc vec array
void PoissonSolve::copy_dat_to_vec(op_dat dat, double *dat_d) {
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, 10, "double", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 1, copy_args);
  memcpy(dat_d, dat->data, dat->set->size * 10 * sizeof(double));
  op_mpi_set_dirtybit(1, copy_args);
}

// Create a PETSc vector for CPUs
void PoissonSolve::create_vec(Vec *v, int size) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, size * mesh->cells->size, PETSC_DECIDE);
}

// Destroy a PETSc vector
void PoissonSolve::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void PoissonSolve::load_vec(Vec *v, op_dat v_dat, int size) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, size, "double", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 1, vec_petsc_args);
  memcpy(v_ptr, (double *)v_dat->data, size * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
  VecRestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void PoissonSolve::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, 10, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 1, vec_petsc_args);
  memcpy((double *)v_dat->data, v_ptr, 10 * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
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

void PoissonSolve::create_shell_mat(Mat *m) {
  MatCreateShell(PETSC_COMM_WORLD, 10 * mesh->cells->size, 10 * mesh->cells->size, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(*m, VECSTANDARD);
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
  global_ind = get_global_start_index(glb_ind->set);
  #endif
  op_arg args[] = {
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 1, args);
  int *data_ptr = (int *)glb_ind->data;
  for(int i = 0; i < mesh->cells->size; i++) {
    data_ptr[i] = global_ind + i;
  }
  op_mpi_set_dirtybit(1, args);
}

void PoissonSolve::setMatrix() {
  if(matCreated) {
    MatDestroy(&Amat);
  }
  MatCreate(PETSC_COMM_WORLD, &Amat);
  matCreated = true;
  MatSetSizes(Amat, 10 * mesh->cells->size, 10 * mesh->cells->size, PETSC_DECIDE, PETSC_DECIDE);

  #ifdef INS_MPI
  MatSetType(Amat, MATMPIAIJ);
  MatMPIAIJSetPreallocation(Amat, 10 * 4, NULL, 0, NULL);
  #else
  MatSetType(Amat, MATSEQAIJ);
  MatSeqAIJSetPreallocation(Amat, 10 * 4, NULL);
  #endif
  MatSetOption(Amat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(op1, -1, OP_ID, 10 * 10, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 2, args);
  double *op1_data = (double *)op1->data;
  int *glb = (int *)glb_ind->data;

  for(int i = 0; i < mesh->cells->size; i++) {
    int global_ind = glb[i];
    // Convert data to row major format
    for(int m = 0; m < 10; m++) {
      for(int n = 0; n < 10; n++) {
        int row = global_ind * 10 + m;
        int col = global_ind * 10 + n;
        double val = op1_data[i * 10 * 10 + m * 10 + n];
        MatSetValues(Amat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }
  }

  op_mpi_set_dirtybit(2, args);

  op_arg edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, 10 * 10, "double", OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, 10 * 10, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->edges, 4, edge_args);

  double *op2L_data = (double *)op2[0]->data;
  double *op2R_data = (double *)op2[1]->data;
  int *glb_l = (int *)glb_indL->data;
  int *glb_r = (int *)glb_indR->data;

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < mesh->edges->size; i++) {
    int leftElement = glb_l[i];
    int rightElement = glb_r[i];

    // Gauss OPf
    // Convert data to row major format
    for(int m = 0; m < 10; m++) {
      for(int n = 0; n < 10; n++) {
        int row = leftElement * 10 + m;
        int col = rightElement * 10 + n;
        double val = op2L_data[i * 10 * 10 + m * 10 + n];
        MatSetValues(Amat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }

    // Convert data to row major format
    for(int m = 0; m < 10; m++) {
      for(int n = 0; n < 10; n++) {
        int row = rightElement * 10 + m;
        int col = leftElement * 10 + n;
        double val = op2R_data[i * 10 * 10 + m * 10 + n];
        MatSetValues(Amat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }
  }

  op_mpi_set_dirtybit(4, edge_args);

  MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY);
}
