#ifndef __INS_BLAS_CALLS_H
#define __INS_BLAS_CALLS_H

#include "op_seq.h"
#include "ins_data.h"
#include "poisson.h"
#include "dg_constants.h"
#include "dg_mesh.h"

extern DGConstants *constants;

void init_gauss_grad_blas(DGMesh *mesh, INSData *data);

void init_gauss_grad_neighbour_blas(DGMesh *mesh, INSData *data);

#endif
