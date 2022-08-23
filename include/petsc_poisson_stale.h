#ifndef __PETSC_POISSON_STALE_H
#define __PETSC_POISSON_STALE_H

#include "op_seq.h"
#include "ins_data.h"
#include "petscvec.h"
#include "petscksp.h"
#include "timing.h"
#include "dg_mesh.h"
#include "ls.h"
#include "poisson_mat.h"

class PetscPoissonStaleSolve {
public:
  PetscPoissonStaleSolve(DGMesh *m, INSData *nsData, LS *s);
  ~PetscPoissonStaleSolve();

  void init();
  bool solve(op_dat b_dat, op_dat x_dat);
  void calc_rhs(const double *u_d, double *rhs_d);
  void precond(const double *in_d, double *out_d);

  double getAverageConvergeIter();

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);

  // OP2 Dats
  op_dat u, rhs, in, out, pre;
  op_dat factor, gFactor, cFactor, mmFactor, h, gDelta;

protected:
  void setMatrix();

  DGMesh *mesh;
  INSData *data;
  LS *ls;

  int dirichlet[3];
  int neumann[3];

  bool massMat;
  double massFactor;
  bool block_jacobi_pre;
  bool pMatInit, vec_created;

  Mat pMat;
  KSP ksp;

  PoissonMat *mat;

  int num_iter_since_refresh;
  int refresh_every;

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

  double *u_data, *rhs_data, *in_data, *out_data, *pre_data;
  double *factor_data, *gFactor_data, *cFactor_data, *mmFactor_data, *h_data;
  double *gDelta_data;
};

class PetscPressureStaleSolve : public PetscPoissonStaleSolve {
public:
  PetscPressureStaleSolve(DGMesh *m, INSData *d, LS *s);

  void setup();
  void set_refresh_num(int rn);
};

#endif
