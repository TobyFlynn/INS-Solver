#ifndef __BLAS_CALLS_H
#define __BLAS_CALLS_H

#include "ins_data.h"

extern double ones[15];
extern double r[15];
extern double s[15];
extern double Dr[15 * 15];
extern double Ds[15 * 15];
extern double Drw[15 * 15];
extern double Dsw[15 * 15];
extern double LIFT[15 * 15];

void init_grid_blas(INSData *nsData);

#endif
