#ifndef __INS_PETSC_UTILS_H
#define __INS_PETSC_UTILS_H

#include "op_seq.h"
#include "petscvec.h"

namespace PETScUtils {
  // General functions
  void create_vec(Vec *v, int local_unknowns);
  void destroy_vec(Vec *v);

  // Not p-adaptive functions

  void dat_to_vec(Vec *v, op_dat v_dat);
  void vec_to_dat(Vec *v, op_dat v_dat);

  void dat_to_ptr(op_dat dat, double *dat_d);
  void ptr_to_dat(op_dat dat, const double *dat_d);

  // p-adaptive functions

  void dat_to_vec(Vec *v, op_dat v_dat, op_dat order);
  void vec_to_dat(Vec *v, op_dat v_dat, op_dat order);

  void dat_to_ptr(op_dat dat, double *dat_d, op_dat order);
  void ptr_to_dat(op_dat dat, const double *dat_d, op_dat order);
}

#endif
