#ifndef __LINEAR_SOLVER_H
#define __LINEAR_SOLVER_H

#include "dg_compiler_defs.h"

#include "op_seq.h"
#include "matrices/poisson_matrix.h"

#include <vector>

class LinearSolver {
public:
  virtual void set_matrix(PoissonMatrix *mat);
  void set_bcs(op_dat bcs);
  void set_nullspace(bool ns);
  virtual bool solve(op_dat rhs, op_dat ans) = 0;
  virtual void init();

protected:
  PoissonMatrix *matrix;
  bool nullspace;
  op_dat bc;
  std::vector<int> iter_counts;
};

#endif
