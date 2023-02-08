#ifndef __PETSC_AMG_H
#define __PETSC_AMG_H

#include "op_seq.h"
#include "linear_solver.h"
#include "petscvec.h"
#include "petscksp.h"

class PETScAMGSolver : public LinearSolver {
public:
  PETScAMGSolver();
  ~PETScAMGSolver();

  bool solve(op_dat rhs, op_dat ans) override;

private:
  KSP ksp;

  bool pMatInit;
  Mat *pMat;
};

#endif
