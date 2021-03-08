#ifndef __BLAS_CALLS_H
#define __BLAS_CALLS_H

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
extern double visMat[15 * 15];
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

void init_cubature_grad_blas(INSData *nsData, CubatureData *cubData);

void init_cubature_blas(INSData *nsData, CubatureData *cubData);

void cubature_op_blas(INSData *nsData, CubatureData *cubData);

void cubature_mm_blas(INSData *nsData, CubatureData *cubData);

void init_gauss_grad_blas(INSData *nsData, GaussData *gaussData);

void init_gauss_blas(INSData *nsData, GaussData *gaussData);

void init_grid_blas(INSData *nsData);

void gauss_op_blas(INSData *nsData, GaussData *gaussData);

void div_blas(INSData *nsData, op_dat u, op_dat v);

void grad_blas(INSData *nsData, op_dat u);

void advection_lift_blas(INSData *nsData, int ind);

void pressure_rhs_blas(INSData *nsData, int ind);

void viscosity_rhs_blas(INSData *nsData);

void poisson_rhs_blas1(INSData *nsData, Poisson *pData);

void poisson_rhs_blas2(INSData *nsData, Poisson *pData);

void poisson_test_rhs_blas(INSData *nsData, op_dat rhs);

void poisson_mass_blas(INSData *nsData, Poisson *pData, double factor);

#endif
