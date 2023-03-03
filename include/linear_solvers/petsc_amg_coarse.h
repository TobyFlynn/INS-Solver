#ifndef __PETSC_AMG_COARSE_H
#define __PETSC_AMG_COARSE_H

#include "op_seq.h"
#include "petsc_amg.h"
#include "petscvec.h"
#include "petscksp.h"
#include "dg_mesh/dg_mesh.h"

class PETScAMGCoarseSolver : public PETScAMGSolver {
public:
  PETScAMGCoarseSolver(DGMesh *m);

  virtual bool solve(op_dat rhs, op_dat ans) override;
};

#endif
