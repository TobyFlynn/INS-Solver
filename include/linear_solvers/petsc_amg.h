#ifndef __PETSC_AMG_H
#define __PETSC_AMG_H

#include "op_seq.h"
#include "linear_solver.h"
#include "petscvec.h"
#include "petscksp.h"
#include "dg_mesh/dg_mesh_2d.h"

class PETScAMGSolver : public LinearSolver {
public:
  PETScAMGSolver(DGMesh2D *m);
  ~PETScAMGSolver();

  bool solve(op_dat rhs, op_dat ans) override;

private:
  DGMesh2D *mesh;
  KSP ksp;

  bool pMatInit;
  Mat *pMat;
};

#endif
