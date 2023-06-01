#ifndef __INS_P_MULTIGRID_H
#define __INS_P_MULTIGRID_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh.h"
#include "petscvec.h"
#include "petscksp.h"
#include "linear_solver.h"
#include "matrices/poisson_coarse_matrix.h"
#include "matrices/poisson_semi_matrix_free.h"
#include "matrices/poisson_matrix_free_diag.h"

#include <vector>

class PMultigridPoissonSolver : public LinearSolver {
public:
  PMultigridPoissonSolver(DGMesh *m);
  ~PMultigridPoissonSolver();

  void init() override;

  virtual void set_matrix(PoissonMatrix *mat) override;
  void set_coarse_matrix(PoissonCoarseMatrix *c_mat);
  bool solve(op_dat rhs, op_dat ans) override;

  void calc_rhs(const DG_FP *u_d, DG_FP *rhs_d);
private:
  void cycle(int order, const int level);
  void smooth(const int iter, const int level);
  void jacobi_smoother(const int level);
  void chebyshev_smoother(const int level);

  DG_FP maxEigenValue();
  void setRandomVector(op_dat vec);
  void setupDirectSolve();

  enum Smoothers {
    JACOBI, CHEBYSHEV
  };

  enum CoarseSolvers {
    PETSC, AMGX, HYPRE
  };

  DGMesh *mesh;
  LinearSolver *coarseSolver;

  PoissonSemiMatrixFree *smfMatrix;
  PoissonMatrixFreeDiag *mfdMatrix;
  PoissonCoarseMatrix *coarseMatrix;
  bool coarseMatCalcRequired;
  bool diagMat;

  Smoothers smoother;
  CoarseSolvers coarseSolver_type;

  // op_dat rk[3], rkQ;

  std::vector<int> orders;
  std::vector<int> pre_it;
  std::vector<int> post_it;
  std::vector<double> eig_vals;
  std::vector<op_dat> u_dat, b_dat, diag_dats, eigen_tmps;
  double eigen_val_saftey_factor;

  int num_levels;

  DG_FP w, max_eig;
};

#endif
