#ifndef __INS_POISSON_H
#define __INS_POISSON_H

#include "ins_data.h"
#include "op_seq.h"
#include "petscvec.h"
#include "petscksp.h"

#include "dg_global_constants.h"
#include "dg_mesh.h"

class PoissonSolve {
public:
  PoissonSolve(DGMesh *m, INSData *d, bool p);
  ~PoissonSolve();

  void init();
  bool solve(op_dat b_dat, op_dat x_dat);
  void calc_rhs(const double *u_d, double *rhs_d);

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);
  double getAverageConvergeIter();

  op_dat u, rhs, h, op1, op2[2], op_bc;
  op_dat factor, gFactor, cFactor, mmFactor;
  op_dat glb_ind, glb_indL, glb_indR, glb_indBC;

protected:
  void set_op();
  void setMatrix();

  DGMesh *mesh;
  INSData *data;
  bool massMat, precondition;
  Mat Amat;
  KSP ksp;

private:
  void create_vec(Vec *v, int size = 15);
  void destroy_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat, int size = 15);
  void store_vec(Vec *v, op_dat v_dat);
  void copy_u(const double *u_d);
  void copy_rhs(double *rhs_d);
  void create_shell_mat(Mat *m);
  void setGlbInd();

  int dirichlet[3], neumann[3];

  int numberIter = 0;
  int solveCount = 0;
  bool matCreated = false;

  op_dat bc_dat;

  Vec b, x;

  double *u_data, *rhs_data, *h_data, *op1_data, *op2_data[2], *op_bc_data;
  double *factor_data, *gFactor_data, *cFactor_data, *mmFactor_data;
  int *glb_ind_data, *glb_indL_data, *glb_indR_data, *glb_indBC_data;
};

class PressureSolve : public PoissonSolve {
public:
  PressureSolve(DGMesh *m, INSData *d, bool p = false) : PoissonSolve(m, d, p) {}

  void setup();
};

class ViscositySolve : public PoissonSolve {
public:
  ViscositySolve(DGMesh *m, INSData *d, bool p = false) : PoissonSolve(m, d, p) {}

  void setup(double mmConst);
};

#endif
