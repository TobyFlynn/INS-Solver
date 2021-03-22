#ifndef __AIRFOIL_SAVESOLUTION_H
#define __AIRFOIL_SAVESOLUTION_H

#include <string>

#include "ins_data.h"

void save_solution(std::string filename, int numPts, int numCells, double *q0,
                   double *q1, double *p, int *cellMap);

void save_solution_cell(std::string filename, int numPts, int numCells, double *q0,
                        double *q1, int *cellMap);

void save_solution_all(std::string filename, INSData *data, int ind);

void save_solution_t(std::string filename, int numPts, int numCells, double *q0, double *q1,
                     double *p, double *pRHS, double *px, double *py, double *utx, double *uty,
                     double *uttx, double *utty, double *visx, double *visy, int *cellMap);

#endif
