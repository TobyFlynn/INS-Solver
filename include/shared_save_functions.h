#ifndef __INS_SHARED_SAVE_FUNCTIONS_H
#define __INS_SHARED_SAVE_FUNCTIONS_H

#include <vector>
#ifdef INS_MPI
#include "pcgnslib.h"
#else
#include "cgnslib.h"
#endif

void get_data_vectors_order_4(std::vector<double> &x_v, std::vector<double> &y_v,
                              std::vector<double> &u_v, std::vector<double> &v_v,
                              std::vector<double> &pr_v, std::vector<double> &vort_v,
                              std::vector<double> &s_v, std::vector<cgsize_t> &cells,
                              double *Ux, double *Uy, double *pr, double *vort,
                              double *x, double *y, double *s, int numCells);

void get_data_vectors_order_3(std::vector<double> &x_v, std::vector<double> &y_v,
                              std::vector<double> &u_v, std::vector<double> &v_v,
                              std::vector<double> &pr_v, std::vector<double> &vort_v,
                              std::vector<double> &s_v, std::vector<cgsize_t> &cells,
                              double *Ux, double *Uy, double *pr, double *vort,
                              double *x, double *y, double *s, int numCells);

void get_data_vectors_order_2(std::vector<double> &x_v, std::vector<double> &y_v,
                              std::vector<double> &u_v, std::vector<double> &v_v,
                              std::vector<double> &pr_v, std::vector<double> &vort_v,
                              std::vector<double> &s_v, std::vector<cgsize_t> &cells,
                              double *Ux, double *Uy, double *pr, double *vort,
                              double *x, double *y, double *s, int numCells);

void get_data_vectors_order_1(std::vector<double> &x_v, std::vector<double> &y_v,
                              std::vector<double> &u_v, std::vector<double> &v_v,
                              std::vector<double> &pr_v, std::vector<double> &vort_v,
                              std::vector<double> &s_v, std::vector<cgsize_t> &cells,
                              double *Ux, double *Uy, double *pr, double *vort,
                              double *x, double *y, double *s, int numCells);

#endif
