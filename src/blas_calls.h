#ifndef __INS_BLAS_CALLS_H
#define __INS_BLAS_CALLS_H

#include "op_seq.h"
#include "ins_data.h"
#include "poisson.h"
#include "constants.h"

extern Constants *constants;

void cubature_op_blas(INSData *nsData, CubatureData *cubData);

void cubature_mm_blas(INSData *nsData, CubatureData *cubData);

void init_gauss_grad_blas(INSData *nsData, GaussData *gaussData);

void init_gauss_blas(INSData *nsData, GaussData *gaussData);

void init_gauss_grad_neighbour_blas(INSData *nsData, GaussData *gaussData);

void init_grid_blas(INSData *nsData);

void gauss_op_blas(INSData *nsData, GaussData *gaussData);

void gauss_opf_blas(INSData *nsData, GaussData *gaussData);

void viscosity_rhs_blas(INSData *nsData, CubatureData *cubatureData);

void poisson_rhs_mass_blas(INSData *data, CubatureData *cubatureData, Poisson_MF *poisson, double factor);

void poisson_mf2_blas(INSData *data, Poisson_MF2 *poisson, CubatureData *cubatureData, bool massMat, double massFactor);

// Assumes matrix is in column major form and both op_dat are defined on the same set
void op2_gemv(bool transpose, int m, int n, double alpha, double *A_ptr, int lda, op_dat x, double beta, op_dat y);

#endif
