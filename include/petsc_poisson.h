#ifndef __PETSC_POISSON_H
#define __PETSC_POISSON_H

#include "op_seq.h"
#include "ins_data.h"
#include "petscvec.h"
#include "petscksp.h"
#include "timing.h"
#include "dg_mesh.h"
#include "ls.h"
#include "poisson_mat.h"

class PetscPoissonSolve {
public:
  PetscPoissonSolve(DGMesh *m, INSData *nsData, LS *s);
  ~PetscPoissonSolve();

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
  void create_shell_mat();
  void set_shell_pc(PC pc);

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
  PetscPressureSolve(DGMesh *m, INSData *d, LS *s);

  void setup();
};

class PetscViscositySolve : public PetscPoissonSolve {
public:
  PetscViscositySolve(DGMesh *m, INSData *d, LS *s);

  void setup(double mmConst);
};

#endif
