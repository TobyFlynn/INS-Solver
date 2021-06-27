#include "pressure_solve.h"
#include "viscosity_solve.h"

#include "op_seq.h"

//  PRESSURE STUFF

// Copy u PETSc vec array to OP2 dat (TODO avoid this copy)
void PressureSolve::copy_u(const double *u_d) {
  op_arg u_copy_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 1, u_copy_args);
  memcpy(u->data, u_d, u->set->size * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, u_copy_args);
}

// Copy rhs OP2 dat to PETSc vec array (TODO avoid this copy)
void PressureSolve::copy_rhs(double *rhs_d) {
  op_arg rhs_copy_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges(data->cells, 1, rhs_copy_args);
  memcpy(rhs_d, rhs->data, rhs->set->size * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, rhs_copy_args);
}

// Create a PETSc vector for CPUs
void PressureSolve::create_vec(Vec *v, int size) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, size * data->cells->size, PETSC_DECIDE);
}

// Destroy a PETSc vector
void PressureSolve::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void PressureSolve::load_vec(Vec *v, op_dat v_dat, int size) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, size, "double", OP_READ)
  };
  op_mpi_halo_exchanges(data->cells, 1, vec_petsc_args);
  memcpy(v_ptr, (double *)v_dat->data, size * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
  VecRestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void PressureSolve::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 1, vec_petsc_args);
  memcpy((double *)v_dat->data, v_ptr, 15 * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
  VecRestoreArrayRead(*v, &v_ptr);
}

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  PressureSolve *poisson;
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

void PressureSolve::create_shell_mat(Mat *m) {
  MatCreateShell(PETSC_COMM_WORLD, 15 * data->cells->size, 15 * data->cells->size, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(*m, VECSTANDARD);
}

//  VISCOSITY STUFF
// Copy u PETSc vec array to OP2 dat (TODO avoid this copy)
void ViscositySolve::copy_u(const double *u_d) {
  op_arg u_copy_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 1, u_copy_args);
  memcpy(u->data, u_d, u->set->size * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, u_copy_args);
}

// Copy rhs OP2 dat to PETSc vec array (TODO avoid this copy)
void ViscositySolve::copy_rhs(double *rhs_d) {
  op_arg rhs_copy_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges(data->cells, 1, rhs_copy_args);
  memcpy(rhs_d, rhs->data, rhs->set->size * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, rhs_copy_args);
}

// Create a PETSc vector for CPUs
void ViscositySolve::create_vec(Vec *v, int size) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, size * data->cells->size, PETSC_DECIDE);
}

// Destroy a PETSc vector
void ViscositySolve::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void ViscositySolve::load_vec(Vec *v, op_dat v_dat, int size) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, size, "double", OP_READ)
  };
  op_mpi_halo_exchanges(data->cells, 1, vec_petsc_args);
  memcpy(v_ptr, (double *)v_dat->data, size * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
  VecRestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void ViscositySolve::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 1, vec_petsc_args);
  memcpy((double *)v_dat->data, v_ptr, 15 * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
  VecRestoreArrayRead(*v, &v_ptr);
}

PetscErrorCode matAMultV(Mat A, Vec x, Vec y) {
  ViscositySolve *poisson;
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

void ViscositySolve::create_shell_mat(Mat *m) {
  MatCreateShell(PETSC_COMM_WORLD, 15 * data->cells->size, 15 * data->cells->size, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMultV);
  MatShellSetVecType(*m, VECSTANDARD);
}
