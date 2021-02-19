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

void init_grid_blas(INSData *nsData);

void div_blas(INSData *nsData, op_dat u, op_dat v);

void grad_blas(INSData *nsData, op_dat u);

void advection_lift_blas(INSData *nsData, int ind);

void pressure_rhs_blas(INSData *nsData, int ind);

void viscosity_rhs_blas(INSData *nsData);

void poisson_rhs_blas1(INSData *nsData, Poisson *pData);

void poisson_rhs_blas2(INSData *nsData, Poisson *pData);

void poisson_test_rhs_blas(INSData *nsData, op_dat rhs);

#endif
