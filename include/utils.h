#ifndef __INS_UTILS_H
#define __INS_UTILS_H

#include "op_seq.h"

double *getOP2Array(op_dat dat);

bool is_point_in_cell(const double x, const double y, const double *cellX, const double *cellY);

void newton_method(const int numPts, double *closest_x, double *closest_y, const double *x, const double *y,
                   int *cell_ind, op_map edge_map, op_dat x_dat, op_dat y_dat, op_dat s_dat, const double h);

#endif
