#include "poisson.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

// Copy u PETSc vec array to OP2 dat (TODO avoid this copy)
void Poisson_MF2::copy_u(const double *u_d) {
  op_arg u_copy_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 1, u_copy_args);
  memcpy(u->data, u_d, u->set->size * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, u_copy_args);
}

// Copy rhs OP2 dat to PETSc vec array (TODO avoid this copy)
void Poisson_MF2::copy_rhs(double *rhs_d) {
  op_arg rhs_copy_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 1, rhs_copy_args);
  memcpy(rhs_d, rhs->data, rhs->set->size * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, rhs_copy_args);
}

// Create a PETSc vector for CPUs
void Poisson::create_vec(Vec *v, int size) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, size * mesh->cells->size, PETSC_DECIDE);
}

// Destroy a PETSc vector
void Poisson::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void Poisson::load_vec(Vec *v, op_dat v_dat, int size) {
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
void Poisson::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 1, vec_petsc_args);
  memcpy((double *)v_dat->data, v_ptr, 15 * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
  VecRestoreArrayRead(*v, &v_ptr);
}

// Create a PETSc matrix for CPUs
void Poisson::create_mat(Mat *m, int row, int col, int prealloc0, int prealloc1) {
  MatCreate(PETSC_COMM_WORLD, m);
  MatSetSizes(*m, row, col, PETSC_DECIDE, PETSC_DECIDE);

  #ifdef INS_MPI
  MatSetType(*m, MATMPIAIJ);
  MatMPIAIJSetPreallocation(*m, prealloc0, NULL, prealloc1, NULL);
  #else
  MatSetType(*m, MATSEQAIJ);
  MatSeqAIJSetPreallocation(*m, prealloc0, NULL);
  #endif
  MatSetOption(*m, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
}

PetscErrorCode matAMult2(Mat A, Vec x, Vec y) {
  timer->startLinearSolveMFMatMult();
  Poisson_MF2 *poisson;
  MatShellGetContext(A, &poisson);
  const double *x_ptr;
  double *y_ptr;
  VecGetArrayRead(x, &x_ptr);
  VecGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecRestoreArrayRead(x, &x_ptr);
  VecRestoreArray(y, &y_ptr);
  timer->endLinearSolveMFMatMult();
  return 0;
}

void Poisson_MF2::create_shell_mat(Mat *m) {
  MatCreateShell(PETSC_COMM_WORLD, 15 * mesh->cells->size, 15 * mesh->cells->size, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMult2);
  MatShellSetVecType(*m, VECSTANDARD);
}

void Poisson_M::setGlbInd() {
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

void Poisson_M::createMassMatrix() {
  create_mat(&pMMat, 15 * mesh->cells->size, 15 * mesh->cells->size, 15);
  pMMatInit = true;

  // Add Cubature OP to mass matrix
  op_arg args[] = {
    op_arg_dat(cData->mm, -1, OP_ID, 15 * 15, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 2, args);
  double *cub_MM = (double *)cData->mm->data;
  int *glb       = (int *)glb_ind->data;

  for(int i = 0; i < mesh->cells->size; i++) {
    int global_ind = glb[i];
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = global_ind * 15 + m;
        int col = global_ind * 15 + n;
        int colInd = n * 15 + m;
        double val = cub_MM[i * 15 * 15 + colInd];
        MatSetValues(pMMat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }
  }

  op_mpi_set_dirtybit(2, args);

  MatAssemblyBegin(pMMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMMat, MAT_FINAL_ASSEMBLY);
}

void Poisson_M::createMatrix() {
  create_mat(&pMat, 15 * mesh->cells->size, 15 * mesh->cells->size, 15 * 4);
  pMatInit = true;
  double tol = 1e-15;

  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 2, args);
  double *op1_data = (double *)op1->data;
  int *glb = (int *)glb_ind->data;

  for(int i = 0; i < mesh->cells->size; i++) {
    int global_ind = glb[i];
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = global_ind * 15 + m;
        int col = global_ind * 15 + n;
        double val = op1_data[i * 15 * 15 + m * 15 + n];
        MatSetValues(pMat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }
  }

  op_mpi_set_dirtybit(2, args);

  op_arg edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, 15 * 15, "double", OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, 15 * 15, "double", OP_READ),
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
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = leftElement * 15 + m;
        int col = rightElement * 15 + n;
        double val = op2L_data[i * 15 * 15 + m * 15 + n];
        MatSetValues(pMat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }

    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = rightElement * 15 + m;
        int col = leftElement * 15 + n;
        double val = op2R_data[i * 15 * 15 + m * 15 + n];
        MatSetValues(pMat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }
  }

  op_mpi_set_dirtybit(4, edge_args);

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
}

void Poisson_M::createBCMatrix() {
  create_mat(&pBCMat, 15 * mesh->cells->size, 21 * mesh->cells->size, 15);
  pBCMatInit = true;
  double tol = 1e-15;

  op_arg args[] = {
    op_arg_dat(op_bc, -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(glb_indBC, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->bedges, 3, args);
  double *op_data = (double *)op_bc->data;
  int *glb = (int *)glb_indBC->data;
  int *edgeNum = (int *)mesh->bedgeNum->data;

  // Create BCs matrix using Gauss data on boundary edges
  for(int i = 0; i < mesh->bedges->size; i++) {
    int global_ind = glb[i];
    for(int j = 0; j < 7 * 15; j++) {
      int col = global_ind * 21 + edgeNum[i] * 7 + (j % 7);
      int row = global_ind * 15 + (j / 7);
      double val = op_data[i * 7 * 15 + j];
      MatSetValues(pBCMat, 1, &row, 1, &col, &val, INSERT_VALUES);
    }
  }

  op_mpi_set_dirtybit(3, args);

  MatAssemblyBegin(pBCMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pBCMat, MAT_FINAL_ASSEMBLY);
}
