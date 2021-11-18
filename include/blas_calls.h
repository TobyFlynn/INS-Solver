#ifndef __INS_BLAS_CALLS_H
#define __INS_BLAS_CALLS_H

#include "op_seq.h"
#include "ins_data.h"
#include "poisson.h"
#include "constants.h"
#include "dg_mesh.h"

extern Constants *constants;

void init_gauss_grad_blas(DGMesh *mesh, GaussData *gaussData);

void init_gauss_grad_neighbour_blas(DGMesh *mesh, GaussData *gaussData);

void poisson_mf2_blas(INSData *data, Poisson_MF2 *poisson, CubatureData *cubatureData, bool massMat, double massFactor);

#endif
