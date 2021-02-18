#ifndef __AIRFOIL_SAVESOLUTION_H
#define __AIRFOIL_SAVESOLUTION_H

#include <string>

void save_solution(std::string filename, int numPts, int numCells, double *q0,
                   double *q1, int *cellMap);

#endif
