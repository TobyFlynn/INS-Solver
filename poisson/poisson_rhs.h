#ifndef __POISSON_RHS_H
#define __POISSON_RHS_H

#include "op_seq.h"
#include "ins_data.h"

extern int FMASK[15];
extern INSData *data;

void poisson_rhs(const double *u, double *rhs);

// void poisson_bc_blas(INSData *nsData);

#endif
