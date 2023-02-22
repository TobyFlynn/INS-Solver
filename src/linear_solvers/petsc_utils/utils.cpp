#include "linear_solvers/petsc_utils.h"

#include "dg_utils.h"

#include "timing.h"

extern Timing *timer;

void get_num_nodes_petsc_utils(const int N, int *Np, int *Nfp) {
  #if DG_DIM == 2
  DGUtils::numNodes2D(N, Np, Nfp);
  #elif DG_DIM == 3
  DGUtils::numNodes3D(N, Np, Nfp);
  #endif
}

// Copy PETSc vec array to OP2 dat
void PETScUtils::copy_vec_to_dat(op_dat dat, const DG_FP *dat_d) {
  timer->startTimer("PETScUtils - copy_vec_to_dat");
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE)
  };
  op_mpi_halo_exchanges(dat->set, 1, copy_args);

  memcpy(dat->data, dat_d, dat->set->size * DG_NP * sizeof(DG_FP));

  op_mpi_set_dirtybit(1, copy_args);
  timer->endTimer("PETScUtils - copy_vec_to_dat");
}

// Copy OP2 dat to PETSc vec array
void PETScUtils::copy_dat_to_vec(op_dat dat, DG_FP *dat_d) {
  timer->startTimer("PETScUtils - copy_dat_to_vec");
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 1, copy_args);

  memcpy(dat_d, dat->data, dat->set->size * DG_NP * sizeof(DG_FP));

  op_mpi_set_dirtybit(1, copy_args);
  timer->endTimer("PETScUtils - copy_dat_to_vec");
}

// Create a PETSc vector for CPUs
void PETScUtils::create_vec(Vec *v, op_set set) {
  timer->startTimer("PETScUtils - create_vec");
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, set->size * DG_NP, PETSC_DECIDE);
  timer->endTimer("PETScUtils - create_vec");
}

// Destroy a PETSc vector
void PETScUtils::destroy_vec(Vec *v) {
  timer->startTimer("PETScUtils - destroy_vec");
  VecDestroy(v);
  timer->endTimer("PETScUtils - destroy_vec");
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void PETScUtils::load_vec(Vec *v, op_dat v_dat) {
  timer->startTimer("PETScUtils - load_vec");
  DG_FP *v_ptr;
  VecGetArray(*v, &v_ptr);

  copy_dat_to_vec(v_dat, v_ptr);

  VecRestoreArray(*v, &v_ptr);
  timer->endTimer("PETScUtils - load_vec");
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void PETScUtils::store_vec(Vec *v, op_dat v_dat) {
  timer->startTimer("PETScUtils - store_vec");
  const DG_FP *v_ptr;
  VecGetArrayRead(*v, &v_ptr);

  copy_vec_to_dat(v_dat, v_ptr);

  VecRestoreArrayRead(*v, &v_ptr);
  timer->endTimer("PETScUtils - store_vec");
}

// P-Adaptive stuff
// Copy PETSc vec array to OP2 dat
void PETScUtils::copy_vec_to_dat_p_adapt(op_dat dat, const DG_FP *dat_d, DGMesh *mesh) {
  timer->startTimer("PETScUtils - copy_vec_to_dat_p_adapt");
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
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
        DG_FP *block_start_dat_c = (DG_FP *)dat->data + block_start * dat->dim;
        memcpy(block_start_dat_c, dat_d + vec_ind, block_count * DG_NP * sizeof(DG_FP));
        vec_ind += DG_NP * block_count;
      }
      block_count = 0;
    }

    DG_FP *v_c = (DG_FP *)dat->data + i * dat->dim;
    int Np, Nfp;
    get_num_nodes_petsc_utils(N, &Np, &Nfp);

    memcpy(v_c, dat_d + vec_ind, Np * sizeof(DG_FP));
    vec_ind += Np;
  }

  if(block_count != 0) {
    DG_FP *block_start_dat_c = (DG_FP *)dat->data + block_start * dat->dim;
    memcpy(block_start_dat_c, dat_d + vec_ind, block_count * DG_NP * sizeof(DG_FP));
    vec_ind += DG_NP * block_count;
  }

  op_mpi_set_dirtybit(2, copy_args);
  timer->endTimer("PETScUtils - copy_vec_to_dat_p_adapt");
}

// Copy OP2 dat to PETSc vec array
void PETScUtils::copy_dat_to_vec_p_adapt(op_dat dat, DG_FP *dat_d, DGMesh *mesh) {
  timer->startTimer("PETScUtils - copy_dat_to_vec_p_adapt");
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
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
        const DG_FP *block_start_dat_c = (DG_FP *)dat->data + block_start * dat->dim;
        memcpy(dat_d + vec_ind, block_start_dat_c, block_count * DG_NP * sizeof(DG_FP));
        vec_ind += DG_NP * block_count;
      }
      block_count = 0;
    }

    const DG_FP *v_c = (DG_FP *)dat->data + i * dat->dim;
    int Np, Nfp;
    get_num_nodes_petsc_utils(N, &Np, &Nfp);

    memcpy(dat_d + vec_ind, v_c, Np * sizeof(DG_FP));
    vec_ind += Np;
  }

  if(block_count != 0) {
    const DG_FP *block_start_dat_c = (DG_FP *)dat->data + block_start * dat->dim;
    memcpy(dat_d + vec_ind, block_start_dat_c, block_count * DG_NP * sizeof(DG_FP));
    vec_ind += DG_NP * block_count;
  }

  op_mpi_set_dirtybit(2, copy_args);
  timer->endTimer("PETScUtils - copy_dat_to_vec_p_adapt");
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void PETScUtils::load_vec_p_adapt(Vec *v, op_dat v_dat, DGMesh *mesh) {
  timer->startTimer("PETScUtils - load_vec_p_adapt");
  DG_FP *v_ptr;
  VecGetArray(*v, &v_ptr);

  copy_dat_to_vec_p_adapt(v_dat, v_ptr, mesh);

  VecRestoreArray(*v, &v_ptr);
  timer->endTimer("PETScUtils - load_vec_p_adapt");
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void PETScUtils::store_vec_p_adapt(Vec *v, op_dat v_dat, DGMesh *mesh) {
  timer->startTimer("PETScUtils - store_vec_p_adapt");
  const DG_FP *v_ptr;
  VecGetArrayRead(*v, &v_ptr);

  copy_vec_to_dat_p_adapt(v_dat, v_ptr, mesh);

  VecRestoreArrayRead(*v, &v_ptr);
  timer->endTimer("PETScUtils - store_vec_p_adapt");
}

void PETScUtils::create_vec_p_adapt(Vec *v, int local_unknowns) {
  timer->startTimer("PETScUtils - create_vec_p_adapt");
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, local_unknowns, PETSC_DECIDE);
  timer->endTimer("PETScUtils - create_vec_p_adapt");
}
