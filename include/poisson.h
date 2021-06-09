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
  virtual bool solve(op_dat b_dat, op_dat x_dat, bool addMass, op_dat factor) = 0;
  virtual void init() = 0;

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
  void create_mat(Mat *m, int row, int col, int prealloc0, int prealloc1 = 0);
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
  Poisson_M(INSData *data, CubatureData *cubData, GaussData *gaussData);
  ~Poisson_M();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);
  bool solve(op_dat b_dat, op_dat x_dat, bool addMass, op_dat factor);
  void init();

  op_dat glb_ind, glb_indL, glb_indR, glb_indBC, op1, op2[2], op_bc;
private:
  void setGlbInd();
  void setOp();
  void setBCOP();
  void createMatrix();
  void createMassMatrix();
  void createBCMatrix();

  Mat pMat, pBCMat, pMMat, op;
  Vec b, bc, rhs, x;
  KSP ksp;

  bool pMatInit = false;
  bool pMMatInit = false;
  bool pBCMatInit = false;

  int *glb_ind_data;
  int *glb_indL_data;
  int *glb_indR_data;
  int *glb_indBC_data;
  double *op1_data;
  double *op2_data[2];
  double *op_bc_data;
};

class Poisson_MF2 : public Poisson {
public:
  Poisson_MF2(INSData *data, CubatureData *cubData, GaussData *gaussData);
  ~Poisson_MF2();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);
  bool solve(op_dat b_dat, op_dat x_dat, bool addMass, op_dat factor);
  void calc_rhs(const double *u_d, double *rhs_d);
  void init();

  op_dat u, rhs, op1, op2[2], op_bc, u_t, rhs_t;

  void setOp();
  void setBCOP();
private:
  void create_shell_mat(Mat *m);
  void copy_u(const double *u_d);
  void copy_rhs(double *rhs_d);

  double *u_data;
  double *rhs_data;
  double *op1_data;
  double *op2_data[2];
  double *op_bc_data;
  double *u_t_data;
  double *rhs_t_data;

  Mat Amat;
  KSP ksp;
  Vec b, x;
};

#endif
