#ifndef __INS_BLAS_CALLS_H
#define __INS_BLAS_CALLS_H

#include "op_seq.h"
#include "ins_data.h"
#include "poisson.h"

extern double ones[15];
extern double r[15];
extern double s[15];
extern double Dr[15 * 15];
extern double Ds[15 * 15];
extern double Drw[15 * 15];
extern double Dsw[15 * 15];
extern double LIFT[15 * 15];
extern double MASS[15 * 15];
extern double cubDr[46 * 15];
extern double cubDs[46 * 15];
extern double cubV[46 * 15];
extern double cubW[46];
extern double cubVDr[46 * 15];
extern double cubVDs[46 * 15];
extern double gFInterp0[7 * 15];
extern double gFInterp1[7 * 15];
extern double gFInterp2[7 * 15];
extern double gF0Dr[7 * 15];
extern double gF0Ds[7 * 15];
extern double gF1Dr[7 * 15];
extern double gF1Ds[7 * 15];
extern double gF2Dr[7 * 15];
extern double gF2Ds[7 * 15];
extern double gInterp[21 * 15];
extern double gFInterp0R[7 * 15];
extern double gFInterp1R[7 * 15];
extern double gFInterp2R[7 * 15];
extern double gF0DrR[7 * 15];
extern double gF0DsR[7 * 15];
extern double gF1DrR[7 * 15];
extern double gF1DsR[7 * 15];
extern double gF2DrR[7 * 15];
extern double gF2DsR[7 * 15];
extern double invMass[15 * 15];

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

void cub_grad_blas(INSData *data, CubatureData *cubatureData, op_dat u);

void cub_grad_blas2(INSData *data, CubatureData *cubatureData, op_dat ux, op_dat uy);

void poisson_rhs_blas1(INSData *data, Poisson_MF *poisson);

#endif
