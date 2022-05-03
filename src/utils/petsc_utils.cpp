#include "petsc_utils.h"

#include "op_seq.h"
#include "dg_utils.h"

/*
 * General functions
 */

// Create a PETSc Vec that will work with GPUs
void PETScUtils::create_vec(Vec *v, int local_unknowns) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, local_unknowns, PETSC_DECIDE);
}

// Destroy a PETSc Vec
void PETScUtils::destroy_vec(Vec *v) {
  VecDestroy(v);
}

/*
 * Not p-adaptive functions
 */

// Copy OP2 dat to PETSc vec
void PETScUtils::dat_to_vec(Vec *v, op_dat v_dat) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);

  PETScUtils::dat_to_ptr(v_dat, v_ptr);

  VecRestoreArray(*v, &v_ptr);
}

// Copy PETSc vec to OP2 dat
void PETScUtils::vec_to_dat(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);

  PETScUtils::ptr_to_dat(v_dat, v_ptr);

  VecRestoreArrayRead(*v, &v_ptr);
}

// Copy OP2 dat to PETSc vec array
void PETScUtils::dat_to_ptr(op_dat dat, double *dat_d) {
  op_arg op2_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 1, op2_args);

  int setSize = dat->set->size;
  memcpy(dat_d, dat->data, setSize * DG_NP * sizeof(double));

  op_mpi_set_dirtybit(1, op2_args);
}

// Copy PETSc vec to OP2 dat
void PETScUtils::ptr_to_dat(op_dat dat, const double *dat_d) {
  op_arg op2_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(dat->set, 1, op2_args);

  int setSize = dat->set->size;
  memcpy(dat->data, dat_d, setSize * DG_NP * sizeof(double));

  op_mpi_set_dirtybit(1, op2_args);
}

/*
 * p-adaptive functions
 */

// Copy OP2 dat to PETSc vec
void PETScUtils::dat_to_vec(Vec *v, op_dat v_dat, op_dat order) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);

  PETScUtils::dat_to_ptr(v_dat, v_ptr, order);

  VecRestoreArray(*v, &v_ptr);
}

// Copy PETSc vec to OP2 dat
void PETScUtils::vec_to_dat(Vec *v, op_dat v_dat, op_dat order) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);

  PETScUtils::ptr_to_dat(v_dat, v_ptr, order);

  VecRestoreArrayRead(*v, &v_ptr);
}

// Copy OP2 dat to PETSc vec array
void PETScUtils::dat_to_ptr(op_dat dat, double *dat_d, op_dat order) {
  op_arg op2_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 2, op2_args);

  int setSize = dat->set->size;
  const int *p = (int *)order->data;

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

  op_mpi_set_dirtybit(2, op2_args);
}

// Copy PETSc vec to OP2 dat
void PETScUtils::ptr_to_dat(op_dat dat, const double *dat_d, op_dat order) {
  op_arg op2_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_WRITE),
    op_arg_dat(order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 2, op2_args);

  int setSize = dat->set->size;
  const int *p = (int *)order->data;

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

  op_mpi_set_dirtybit(2, op2_args);
}
