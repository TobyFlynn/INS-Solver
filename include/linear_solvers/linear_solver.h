#ifndef __LINEAR_SOLVER_H
#define __LINEAR_SOLVER_H

#include "dg_compiler_defs.h"

#include "op_seq.h"
#include "matrices/poisson_matrix.h"

#if DG_FP == double
#define LIN_SOL_TOL 1e-10
#else
#define LIN_SOL_TOL 1e-6
#endif

class LinearSolver {
public:
  void set_matrix(PoissonMatrix *mat);
  void set_bcs(op_dat bcs);
  void set_nullspace(bool ns);
  virtual bool solve(op_dat rhs, op_dat ans) = 0;

protected:
  PoissonMatrix *matrix;
  bool nullspace;
  op_dat bc;
};

#endif
