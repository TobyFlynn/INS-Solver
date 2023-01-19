#ifndef __INS_P_MULTIGRID_H
#define __INS_P_MULTIGRID_H

#include "op_seq.h"
#include "dg_mesh.h"
#include "poisson_mat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "poisson_mat.h"

class PMultigrid {
public:
  PMultigrid(DGMesh *m);
  ~PMultigrid();

  void init();
  bool solve(op_dat b, op_dat x);
  bool sub_solve(op_dat b, op_dat u);

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);

  void calc_rhs(const double *u_d, double *rhs_d);
  void set_rho(op_dat rho);

private:
  bool cycle(int p);

  double maxEigenValue();
  void setRandomVector(op_dat vec);

  void copy_vec_to_dat(op_dat dat, const double *dat_d);
  void copy_dat_to_vec(op_dat dat, double *dat_d);
  void create_shell_mat(Mat *mat);
  void create_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat);
  void store_vec(Vec *v, op_dat v_dat);
  void setMatrix();

  Mat pMat;
  bool pMatInit, vec_created;

  KSP ksp;
  Vec b, x;

  DGMesh *mesh;
  int dirichlet[3];
  int neumann[3];

  op_dat bc_dat;

  op_dat tmp_dat[DG_ORDER], u_dat[DG_ORDER], b_dat[DG_ORDER], u_rhs, rhs_rhs;
  op_dat fact[DG_ORDER];

  PoissonMat *pMatrix;
};

#endif
