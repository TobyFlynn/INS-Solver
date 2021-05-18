#ifndef __INS_SAVESOLUTION_H
#define __INS_SAVESOLUTION_H

#include <string>

#include "ins_data.h"

void save_solution_init(std::string filename, INSData *data);

void save_solution_iter(std::string filename, INSData *data, int ind, int iter);

void save_solution_finalise(std::string filename, INSData *data, int numIter, double dt);

void save_solution(std::string filename, INSData *data, int ind, double finalTime = 0.0, double nu = 0.0);

#endif
