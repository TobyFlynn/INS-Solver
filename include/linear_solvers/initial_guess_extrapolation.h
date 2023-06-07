#ifndef __INITIAL_GUESS_EXTRAPOLATION_H
#define __INITIAL_GUESS_EXTRAPOLATION_H

#include <vector>

#include "dg_dat_pool.h"
#include "dg_mesh/dg_mesh.h"

#define EXTRAPOLATE_HISTORY_SIZE 8

void initial_guess_extrapolation(DGMesh *mesh, std::vector<std::pair<DG_FP,DGTempDat>> &history,
                                 op_dat init_guess, const DG_FP t_n1);

#endif
