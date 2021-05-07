#ifndef __INS_OPERATORS_H
#define __INS_OPERATORS_H

#include "op_seq.h"
#include "ins_data.h"

void div(INSData *data, op_dat u, op_dat v, op_dat res);

void curl(INSData *data, op_dat u, op_dat v, op_dat res);

void grad(INSData *data, op_dat u, op_dat ux, op_dat uy);

void cub_grad(INSData *data, CubatureData *cubData, op_dat u, op_dat ux, op_dat uy, bool weak = false);

void cub_div(INSData *data, CubatureData *cubData, op_dat u, op_dat v, op_dat res, bool weak = false);

#endif