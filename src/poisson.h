#ifndef __POISSON_H
#define __POISSON_H

#include "op_seq.h"
#include "ins_data.h"
#include "petscvec.h"
#include "petscksp.h"
#include "timing.h"

extern double gFInterp0_g[7 * 15];
extern double gFInterp1_g[7 * 15];
extern double gFInterp2_g[7 * 15];
extern double gaussW_g[7];
extern Timing *timer;

class Poisson {
public:
  Poisson(INSData *nsData, CubatureData *cubData, GaussData *gaussData);
  ~Poisson();

  virtual bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0) = 0;

  double getAverageConvergeIter();

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);

  // OP2 Dats
  op_dat bc_dat;
protected:
  void create_vec(Vec *v, int size = 15);
  void destroy_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat, int size = 15);
  void store_vec(Vec *v, op_dat v_dat);
  void create_mat(Mat *m, int row, int col, int prealloc);
  INSData *data;
  CubatureData *cData;
  GaussData *gData;

  int *dirichlet;
  int *neumann;

  bool massMat;
  double massFactor;

  int numberIter = 0;
  int solveCount = 0;
};

class Poisson_M : public Poisson {
public:
  Poisson_M(INSData *data, CubatureData *cubData, GaussData *gaussData) : Poisson(data, cubData, gaussData) {}
  ~Poisson_M();

  void createMatrix();
  void createMassMatrix();
  void createBCMatrix();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);

private:
  Mat pMat, pBCMat, pMMat;

  bool pMatInit = false;
  bool pMMatInit = false;
  bool pBCMatInit = false;
};

class Poisson_MF : public Poisson {
public:
  Poisson_MF(INSData *data, CubatureData *cubData, GaussData *gaussData);
  ~Poisson_MF();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);
  void calc_rhs(const double *u_d, double *rhs_d);

  op_dat u, rhs, tau, gU, uNumFlux, uFluxX, uFluxY, dudx, dudy, qx, qy, gqx, gqy;
  op_dat qxNumFlux, qyNumFlux, qFlux, gradx, grady;

private:
  void create_shell_mat(Mat *m);
  void copy_u(const double *u_d);
  void copy_rhs(double *rhs_d);
  void applyBCs(op_dat b_dat);

  double *u_data;
  double *rhs_data;
  double *tau_data;
  double *gU_data;
  double *uNumFlux_data;
  double *uFluxX_data;
  double *uFluxY_data;
  double *dudx_data;
  double *dudy_data;
  double *qx_data;
  double *qy_data;
  double *gqx_data;
  double *gqy_data;
  double *qxNumFlux_data;
  double *qyNumFlux_data;
  double *qFlux_data;
  double *gradx_data;
  double *grady_data;
};

#endif
