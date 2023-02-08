#include "linear_solvers/petsc_utils.h"

#include "dg_utils.h"

// Copy PETSc vec array to OP2 dat
void PETScUtils::copy_vec_to_dat(op_dat dat, const double *dat_d) {
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(dat->set, 1, copy_args);

  memcpy(dat->data, dat_d, dat->set->size * DG_NP * sizeof(double));

  op_mpi_set_dirtybit(1, copy_args);
}

// Copy OP2 dat to PETSc vec array
void PETScUtils::copy_dat_to_vec(op_dat dat, double *dat_d) {
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 1, copy_args);

  memcpy(dat_d, dat->data, dat->set->size * DG_NP * sizeof(double));

  op_mpi_set_dirtybit(1, copy_args);
}

// Create a PETSc vector for CPUs
void PETScUtils::create_vec(Vec *v, op_set set) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, set->size * DG_NP, PETSC_DECIDE);
}

// Destroy a PETSc vector
void PETScUtils::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void PETScUtils::load_vec(Vec *v, op_dat v_dat) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);

  copy_dat_to_vec(v_dat, v_ptr);

  VecRestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void PETScUtils::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);

  copy_vec_to_dat(v_dat, v_ptr);

  VecRestoreArrayRead(*v, &v_ptr);
}

// P-Adaptive stuff
// Copy PETSc vec array to OP2 dat
void PETScUtils::copy_vec_to_dat_p_adapt(op_dat dat, const double *dat_d, DGMesh2D *mesh) {
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
    DGUtils::numNodes2D(N, &Np, &Nfp);

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
void PETScUtils::copy_dat_to_vec_p_adapt(op_dat dat, double *dat_d, DGMesh2D *mesh) {
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
    DGUtils::numNodes2D(N, &Np, &Nfp);

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

// Load a PETSc vector with values from an OP2 dat for CPUs
void PETScUtils::load_vec_p_adapt(Vec *v, op_dat v_dat, DGMesh2D *mesh) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);

  copy_dat_to_vec_p_adapt(v_dat, v_ptr, mesh);

  VecRestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void PETScUtils::store_vec_p_adapt(Vec *v, op_dat v_dat, DGMesh2D *mesh) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);

  copy_vec_to_dat_p_adapt(v_dat, v_ptr, mesh);

  VecRestoreArrayRead(*v, &v_ptr);
}

void PETScUtils::create_vec_p_adapt(Vec *v, int local_unknowns) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, local_unknowns, PETSC_DECIDE);
}
