#ifndef __POISSON_H
#define __POISSON_H

#include "op_seq.h"
#include "ins_data.h"
#include "petscvec.h"
#include "petscksp.h"
#include "timing.h"

#include "dg_global_constants.h"
#include "dg_mesh.h"

extern Timing *timer;

class Poisson {
public:
  Poisson(DGMesh *m, INSData *d);
  ~Poisson();

  virtual bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0) = 0;
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
  DGMesh *mesh;
  INSData *data;

  int dirichlet[3];
  int neumann[3];

  bool massMat;
  double massFactor;
  bool scalarFactor;
  op_dat massFactorDat;

  int numberIter = 0;
  int solveCount = 0;
};

class Poisson_M : public Poisson {
public:
  Poisson_M(DGMesh *m, INSData *d);
  ~Poisson_M();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);
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

  int *glb_ind_data, *glb_indL_data, *glb_indR_data, *glb_indBC_data;
  double *op1_data, *op2_data[2], *op_bc_data;
};

class Poisson_MF : public Poisson {
public:
  Poisson_MF(DGMesh *m, INSData *d);
  ~Poisson_MF();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);
  void calc_rhs(const double *u_d, double *rhs_d);
  void init();

  op_dat u, rhs, gU, gDudx, gDudy, fluxX, fluxY, flux, dudx, dudy, qx, qy, tmp_u;

private:
  void create_shell_mat(Mat *m);
  void copy_u(const double *u_d);
  void copy_rhs(double *rhs_d);

  void apply_bc(op_dat b);

  double *u_data, *rhs_data, *gU_data, *gDudx_data, *gDudy_data;
  double *fluxX_data, *fluxY_data, *flux_data, *dudx_data, *dudy_data;
  double *qx_data, *qy_data, *tmp_u_data;

  Mat Amat;
  KSP ksp;
  Vec b, x;
};

class Poisson_MF2 : public Poisson {
public:
  Poisson_MF2(DGMesh *m, INSData *d);
  ~Poisson_MF2();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);
  void calc_rhs(const double *u_d, double *rhs_d);
  void init();

  op_dat u, rhs, op1, op2[2], op_bc, u_t, rhs_t;

  void setOp();
  void setBCOP();
private:
  void create_shell_mat(Mat *m);
  void copy_u(const double *u_d);
  void copy_rhs(double *rhs_d);

  double *u_data, *rhs_data, *op1_data, *op2_data[2], *op_bc_data;
  double *u_t_data, *rhs_t_data;

  Mat Amat;
  KSP ksp;
  Vec b, x;
};

#endif
