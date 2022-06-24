#ifndef __INS_UTILS_H
#define __INS_UTILS_H

#include "op_seq.h"

double *getOP2Array(op_dat dat);

void xy_to_rs(const double x, const double y, double &r, double &s, const double *cellX, const double *cellY);

void rs_to_xy(const double r, const double s, double &x, double &y, const double *cellX, const double *cellY);

double eval_at_pt(const double r, const double s, const double *modal);

void eval_grad_at_pt(const double r, const double s, const double *modal,
                     double &dr, double &ds);

void eval_hessian_at_pt(const double r, const double s, const double *modal,
                        double &dr2, double &drs, double &ds2);

bool is_point_in_cell(const double x, const double y, const double *cellX, const double *cellY);

void newton_method(const int numPts, double *s, double *closest_x,
                   double *closest_y, int *cell_ind, const double *x,
                   const double *y, const double *s_modal, const double *cell_x,
                   const double *cell_y, const double *rx, const double *sx,
                   const double *ry, const double *sy, const double h);

#endif
