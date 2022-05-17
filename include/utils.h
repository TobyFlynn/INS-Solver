#ifndef __INS_UTILS_H
#define __INS_UTILS_H

#include "op_seq.h"

double *getOP2Array(op_dat dat);

double xy_to_rs(const double x, const double y, double &r, double &s, const double *cellX, const double *cellY);

double rs_to_xy(const double r, const double s, double &x, double &y, const double *cellX, const double *cellY);

double eval_at_pt(const double r, const double s, const double *modal);

void eval_grad_at_pt(const double r, const double s, const double *modal,
                     double &dr, double &ds);

void eval_hessian_at_pt(const double r, const double s, const double *modal,
                        double &dr2, double &drs, double &ds2);

bool is_point_in_cell(const double x, const double y, const double *cellX, const double *cellY);

void newton_method(const double node_x, const double node_y,
                   double &closest_pt_x, double &closest_pt_y,
                   const double *s_modal, const double *dsdr_modal,
                   const double *dsds_modal, const double *dsdr2_modal,
                   const double *dsdrs_modal, const double *dsds2_modal,
                   const double *cellX, const double *cellY);

#endif
