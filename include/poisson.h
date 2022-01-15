#ifndef __POISSON_H
#define __POISSON_H

#include "op_seq.h"
#include "ins_data.h"
#include "petscvec.h"
#include "petscksp.h"
#include "timing.h"
#include "dg_mesh.h"
#include "ls.h"

extern Timing *timer;

class PoissonSolve {
public:
  PoissonSolve(DGMesh *m, INSData *nsData, LS *s);
  ~PoissonSolve();

  void init();
  bool solve(op_dat b_dat, op_dat x_dat);
  void calc_rhs(const double *u_d, double *rhs_d);
  void precond(const double *in_d, double *out_d);

  double getAverageConvergeIter();

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);

  // OP2 Dats
  op_dat op1, op2[2], op_bc;
  op_dat glb_ind, glb_indL, glb_indR, glb_indBC;
  op_dat u, rhs, in, out, pre;
  op_dat factor, gFactor, cFactor, mmFactor, h, gDelta;
  op_dat orderL, orderR, orderBC;

  int unknowns;

protected:
  void set_op();
  void setMatrix();
  void create_shell_mat(Mat *m);
  void set_shell_pc(PC pc);
  int get_local_unknowns();
  void update_glb_ind();

  DGMesh *mesh;
  INSData *data;
  LS *ls;

  int dirichlet[3];
  int neumann[3];

  bool massMat;
  double massFactor;
  bool block_jacobi_pre;
  bool pMatInit;

  Mat pMat;
  KSP ksp;

private:
  void create_vec(Vec *v);
  void destroy_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat);
  void store_vec(Vec *v, op_dat v_dat);
  void copy_vec_to_dat(op_dat dat, const double *dat_d);
  void copy_dat_to_vec(op_dat dat, double *dat_d);

  void setGlbInd();

  void calc_cub_sub_mat();
  void calc_gauss_sub_mat();
  void calc_mm_mat();

  op_dat bc_dat;
  Vec b, x;

  int numberIter, solveCount;

  double *op1_data, *op2_data[2], *op_bc_data;
  int *glb_ind_data, *glb_indL_data, *glb_indR_data, *glb_indBC_data;
  double *u_data, *rhs_data, *in_data, *out_data, *pre_data;
  double *factor_data, *gFactor_data, *cFactor_data, *mmFactor_data, *h_data;
  double *gDelta_data;
  int *orderL_data, *orderR_data, *orderBC_data;
};

class PressureSolve : public PoissonSolve {
public:
  PressureSolve(DGMesh *m, INSData *d, LS *s);

  void setup();
};

class ViscositySolve : public PoissonSolve {
public:
  ViscositySolve(DGMesh *m, INSData *d, LS *s);

  void setup(double mmConst);
};

#endif
