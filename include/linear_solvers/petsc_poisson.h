#ifndef __PETSC_POISSON_H
#define __PETSC_POISSON_H

#include "op_seq.h"
#include "petscvec.h"
#include "petscksp.h"
#include "timing.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "matrices/2d/factor_poisson_matrix_2d.h"

class PetscPoissonSolve {
public:
  PetscPoissonSolve(DGMesh2D *m);
  ~PetscPoissonSolve();

  void init();
  bool solve(op_dat b_dat, op_dat x_dat);
  void calc_rhs(const double *u_d, double *rhs_d);
  void precond(const double *in_d, double *out_d);

  double getAverageConvergeIter();

  void setBCValues(op_dat bc);

  // OP2 Dats
  op_dat u, rhs, in, out, pre;
  op_dat factor, gFactor, cFactor, mmFactor, h, gDelta;

protected:
  void setMatrix();
  void create_shell_mat();
  void set_shell_pc(PC pc);

  DGMesh2D *mesh;

  int dirichlet[3];
  int neumann[3];

  bool massMat;
  double massFactor;
  bool block_jacobi_pre;
  bool pMatInit, vec_created;

  Mat pMat;
  KSP ksp;

  FactorPoissonMatrix2D *mat;

private:
  void create_vec(Vec *v);
  void destroy_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat);
  void store_vec(Vec *v, op_dat v_dat);
  void copy_vec_to_dat(op_dat dat, const double *dat_d);
  void copy_dat_to_vec(op_dat dat, double *dat_d);

  op_dat bc_dat;
  Vec b, x;

  int numberIter, solveCount;
  int prev_unknowns;
};

class PetscPressureSolve : public PetscPoissonSolve {
public:
  PetscPressureSolve(DGMesh2D *m);

  void setup(op_dat rho);
};

class PetscViscositySolve : public PetscPoissonSolve {
public:
  PetscViscositySolve(DGMesh2D *m);

  void setup(double mmConst, op_dat rho, op_dat mu);
private:
  void calc_precond_mat();
};

#endif
