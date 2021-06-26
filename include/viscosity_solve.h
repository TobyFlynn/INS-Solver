#ifndef __INS_VISCOSITY_SOLVE_H
#define __INS_VISCOSITY_SOLVE_H

#include "ins_data.h"
#include "op_seq.h"
#include "petscvec.h"
#include "petscksp.h"

class ViscositySolve {
public:
  ViscositySolve(INSData *d, CubatureData *c, GaussData *g);
  ~ViscositySolve();

  void init();
  bool solve(op_dat b_dat, op_dat x_dat, double factor);
  void calc_rhs(const double *u_d, double *rhs_d);

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);
  double getAverageConvergeIter();

  op_dat u, rhs, h, cMu, gMu, cRho, gRho;

private:
  void create_vec(Vec *v, int size = 15);
  void destroy_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat, int size = 15);
  void store_vec(Vec *v, op_dat v_dat);
  void copy_u(const double *u_d);
  void copy_rhs(double *rhs_d);
  void create_shell_mat(Mat *m);

  INSData *data;
  CubatureData *cData;
  GaussData *gData;

  int *dirichlet;
  int *neumann;

  int numberIter = 0;
  int solveCount = 0;
  double massFactor;

  op_dat bc_dat;

  Mat Amat;
  KSP ksp;
  Vec b, x;

  double *u_data, *rhs_data, *h_data, *cMu_data, *gMu_data, *cRho_data, *gRho_data;
};

#endif
