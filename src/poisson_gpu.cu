#include "poisson.h"

#include "op_seq.h"

// Copy u PETSc vec array to OP2 dat (TODO avoid this copy)
void PoissonSolve::copy_u(const double *u_d) {
  op_arg u_copy_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 1, u_copy_args);
  cudaMemcpy(u->data_d, u_d, u->set->size * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, u_copy_args);
}

// Copy rhs OP2 dat to PETSc vec array (TODO avoid this copy)
void PoissonSolve::copy_rhs(double *rhs_d) {
  op_arg rhs_copy_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 1, rhs_copy_args);
  cudaMemcpy(rhs_d, rhs->data_d, rhs->set->size * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, rhs_copy_args);
}

// Create a PETSc vector for GPUs
void PoissonSolve::create_vec(Vec *v, int size) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECCUDA);
  VecSetSizes(*v, size * mesh->cells->size, PETSC_DECIDE);
}

// Destroy a PETSc vector
void PoissonSolve::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for GPUs
void PoissonSolve::load_vec(Vec *v, op_dat v_dat, int size) {
  double *v_ptr;
  VecCUDAGetArray(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, size, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 1, vec_petsc_args);
  cudaMemcpy(v_ptr, (double *)v_dat->data_d, size * v_dat->set->size * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, vec_petsc_args);
  VecCUDARestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for GPUs
void PoissonSolve::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecCUDAGetArrayRead(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 1, vec_petsc_args);
  cudaMemcpy((double *)v_dat->data_d, v_ptr, 15 * v_dat->set->size * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, vec_petsc_args);
  VecCUDARestoreArrayRead(*v, &v_ptr);
}

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

void PoissonSolve::create_shell_mat(Mat *m) {
  MatCreateShell(PETSC_COMM_WORLD, 15 * mesh->cells->size, 15 * mesh->cells->size, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(*m, VECCUDA);
}
