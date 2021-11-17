#ifndef __POISSON_H
#define __POISSON_H

#include "op_seq.h"
#include "ins_data.h"
#include "petscvec.h"
#include "petscksp.h"
#include "timing.h"

#include "dg_mesh.h"

extern Timing *timer;

class Poisson {
public:
  Poisson(INSData *insData, DGMesh *m);
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
  INSData *data;
  DGMesh *mesh;

  int *dirichlet;
  int *neumann;

  bool massMat;
  double massFactor;

  int numberIter = 0;
  int solveCount = 0;
};

class Poisson_MF2 : public Poisson {
public:
  Poisson_MF2(INSData *insData, DGMesh *m);
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
