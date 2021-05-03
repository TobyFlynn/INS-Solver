#ifndef __INS_BLAS_CALLS_H
#define __INS_BLAS_CALLS_H

#include "op_seq.h"
#include "ins_data.h"
#include "poisson.h"
#include "constants.h"

extern Constants *constants;

void init_cubature_grad_blas(INSData *nsData, CubatureData *cubData);

void init_cubature_blas(INSData *nsData, CubatureData *cubData);

void cubature_op_blas(INSData *nsData, CubatureData *cubData);

void cubature_mm_blas(INSData *nsData, CubatureData *cubData);

void init_gauss_coords_blas(INSData *nsData, GaussData *gaussData);

void init_gauss_grad_blas(INSData *nsData, GaussData *gaussData);

void init_gauss_blas(INSData *nsData, GaussData *gaussData);

void init_gauss_grad_neighbour_blas(INSData *nsData, GaussData *gaussData);

void init_grid_blas(INSData *nsData);

void gauss_op_blas(INSData *nsData, GaussData *gaussData);

void gauss_opf_blas(INSData *nsData, GaussData *gaussData);

void div_blas(INSData *nsData, op_dat u, op_dat v);

void grad_blas(INSData *nsData, op_dat u);

void advection_lift_blas(INSData *nsData, int ind);

void pressure_rhs_blas(INSData *nsData, int ind);

void viscosity_rhs_blas(INSData *nsData, CubatureData *cubatureData);

void gauss_interp_blas(INSData *data, op_dat input, op_dat output);

void cub_grad_w_blas(INSData *data, CubatureData *cubatureData, op_dat u);

void cub_grad_w_blas2(INSData *data, CubatureData *cubatureData, op_dat ux, op_dat uy);

void cub_grad_blas(INSData *data, CubatureData *cubatureData, op_dat u);

void cub_grad_blas2(INSData *data, CubatureData *cubatureData, op_dat ux, op_dat uy);

void poisson_rhs_blas1(INSData *data, Poisson_MF *poisson);

void poisson_rhs_blas2(INSData *data, Poisson_MF *poisson);

void cub_div_w_blas(INSData *data, CubatureData *cubatureData, op_dat u, op_dat v);

void cub_div_w_blas2(INSData *data, CubatureData *cubatureData, op_dat res);

void cub_div_blas(INSData *data, CubatureData *cubatureData, op_dat u, op_dat v);

void cub_div_blas2(INSData *data, CubatureData *cubatureData, op_dat res);

void poisson_rhs_mass_blas(INSData *data, CubatureData *cubatureData, Poisson_MF *poisson, double factor);

void poisson_test_rhs_blas(INSData *nsData, op_dat rhs);

void poisson_bc_blas(INSData *data, Poisson_MF *poisson);

void poisson_bc_blas2(INSData *data, Poisson_MF *poisson);

void poisson_mf2_blas(INSData *data, Poisson_MF2 *poisson, CubatureData *cubatureData, bool massMat, double massFactor);

#endif
