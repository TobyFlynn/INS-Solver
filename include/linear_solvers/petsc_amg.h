#ifndef __PETSC_AMG_H
#define __PETSC_AMG_H

#include "op_seq.h"
#include "linear_solver.h"
#include "petscvec.h"
#include "petscksp.h"
#include "dg_mesh/dg_mesh.h"

class PETScAMGSolver : public LinearSolver {
public:
  PETScAMGSolver(DGMesh *m);
  ~PETScAMGSolver();

  bool solve(op_dat rhs, op_dat ans) override;
  void set_tol(const double tol);

private:
  DGMesh *mesh;
  KSP ksp;

  bool pMatInit;
  Mat *pMat;
};

#endif
