#ifndef __INS_OPERATORS_H
#define __INS_OPERATORS_H

#include "op_seq.h"
#include "ins_data.h"

void div(INSData *data, op_dat u, op_dat v, op_dat res);

void curl(INSData *data, op_dat u, op_dat v, op_dat res);

void grad(INSData *data, op_dat u, op_dat ux, op_dat uy);

void cub_grad(INSData *data, CubatureData *cData, op_dat u, op_dat ux, op_dat uy);

void cub_grad_weak(INSData *data, CubatureData *cData, op_dat u, op_dat ux, op_dat uy);

void cub_div_weak(INSData *data, CubatureData *cData, op_dat u, op_dat v, op_dat res);

void inv_mass(INSData *data, op_dat u);

#endif
