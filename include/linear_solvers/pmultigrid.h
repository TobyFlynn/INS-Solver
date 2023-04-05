#ifndef __INS_P_MULTIGRID_H
#define __INS_P_MULTIGRID_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh.h"
#include "petscvec.h"
#include "petscksp.h"
#include "linear_solver.h"
#include "petsc_amg_coarse.h"
#include "matrices/poisson_coarse_matrix.h"
#include "matrices/poisson_semi_matrix_free.h"
#include "matrices/poisson_matrix_free_diag.h"

#include <vector>

class PMultigridPoissonSolver : public LinearSolver {
public:
  PMultigridPoissonSolver(DGMesh *m);
  ~PMultigridPoissonSolver();

  virtual void set_matrix(PoissonMatrix *mat) override;
  void set_coarse_matrix(PoissonCoarseMatrix *c_mat);
  bool solve(op_dat rhs, op_dat ans) override;

  void calc_rhs(const DG_FP *u_d, DG_FP *rhs_d);
private:
  void cycle(int order, const int level);
  void smoother(const int order, const int level);

  DG_FP maxEigenValue();
  void setRandomVector(op_dat vec);
  void setupDirectSolve();

  DGMesh *mesh;
  PETScAMGCoarseSolver *coarseSolver;

  PoissonSemiMatrixFree *smfMatrix;
  PoissonMatrixFreeDiag *mfdMatrix;
  PoissonCoarseMatrix *coarseMatrix;
  bool coarseMatCalcRequired;
  bool diagMat;

  op_dat eg_tmp_0, eg_tmp_1, rk[3], rkQ;

  std::vector<int> orders;
  std::vector<int> pre_it;
  std::vector<int> post_it;
  std::vector<op_dat> tmp_dat, u_dat, b_dat;
  int num_eigen_val_iter;
  double coarse_solve_tol;

  int num_levels;

  DG_FP w;
};

#endif
