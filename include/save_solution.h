#ifndef __INS_SAVESOLUTION_H
#define __INS_SAVESOLUTION_H

#include <string>

#include "dg_mesh.h"
#include "ins_data.h"
#include "ls.h"

void save_solution_init(std::string filename, DGMesh *mesh, INSData *data, LS *ls);

void save_solution_iter(std::string filename, DGMesh *mesh, INSData *data, int ind, LS *ls, int iter);

void save_solution_finalise(std::string filename, int numIter, double dt);

void save_solution(std::string filename, DGMesh *mesh, INSData *data, int ind, LS *ls, double finalTime = 0.0, double nu = 0.0);

#endif
