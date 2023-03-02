#ifndef __PETSC_UTILS_H
#define __PETSC_UTILS_H

#include "op_seq.h"

#include "petscvec.h"

#include "dg_mesh/dg_mesh.h"

namespace PETScUtils {

  // Copy PETSc vec array to OP2 dat
  void copy_vec_to_dat(op_dat dat, const DG_FP *dat_d);

  // Copy OP2 dat to PETSc vec array
  void copy_dat_to_vec(op_dat dat, DG_FP *dat_d);

  // Create a PETSc vector for CPUs
  void create_vec(Vec *v, op_set set);

  // Destroy a PETSc vector
  void destroy_vec(Vec *v);

  // Load a PETSc vector with values from an OP2 dat for CPUs
  void load_vec(Vec *v, op_dat v_dat);

  // Load an OP2 dat with the values from a PETSc vector for CPUs
  void store_vec(Vec *v, op_dat v_dat);

  // P-Adaptive stuff
  // Copy PETSc vec array to OP2 dat
  void copy_vec_to_dat_p_adapt(op_dat dat, const DG_FP *dat_d, DGMesh *mesh);

  // Copy OP2 dat to PETSc vec array
  void copy_dat_to_vec_p_adapt(op_dat dat, DG_FP *dat_d, DGMesh *mesh);

  // Load a PETSc vector with values from an OP2 dat for CPUs
  void load_vec_p_adapt(Vec *v, op_dat v_dat, DGMesh *mesh);

  // Load an OP2 dat with the values from a PETSc vector for CPUs
  void store_vec_p_adapt(Vec *v, op_dat v_dat, DGMesh *mesh);

  void create_vec_p_adapt(Vec *v, int local_unknowns);

  // Coarse stuff
  // Copy PETSc vec array to OP2 dat
  void copy_vec_to_dat_coarse(op_dat dat, const DG_FP *dat_d);

  // Copy OP2 dat to PETSc vec array
  void copy_dat_to_vec_coarse(op_dat dat, DG_FP *dat_d);

  // Load a PETSc vector with values from an OP2 dat for CPUs
  void load_vec_coarse(Vec *v, op_dat v_dat);

  // Load an OP2 dat with the values from a PETSc vector for CPUs
  void store_vec_coarse(Vec *v, op_dat v_dat);
}

#endif
