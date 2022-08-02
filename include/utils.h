#ifndef __INS_UTILS_H
#define __INS_UTILS_H

#include "op_seq.h"

double *getOP2Array(op_dat dat);

bool is_point_in_cell(const double x, const double y, const double *cellX, const double *cellY);

#endif
