#ifndef __OPERATORS_H
#define __OPERATORS_H

#include "op_seq.h"
#include "ins_data.h"

void div(INSData *data, op_dat u, op_dat v, op_dat res);

void curl(INSData *data, op_dat u, op_dat v, op_dat res);

void grad(INSData *data, op_dat u, op_dat ux, op_dat uy);

#endif
