#ifndef __INS_P_MULTIGRID_H
#define __INS_P_MULTIGRID_H

#include "op_seq.h"
#include "dg_mesh.h"
#include "ins_data.h"
#include "poisson_mat.h"
#include "petscvec.h"
#include "petscksp.h"

class PMultigrid {
public:
  PMultigrid(DGMesh *m, INSData *insData);
  ~PMultigrid();

  void init();
  bool solve(op_dat b, op_dat x);
  void sub_solve(op_dat b, op_dat u);

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);

  void calc_rhs(const double *u_d, double *rhs_d);

  int unknowns;

  op_dat factor, mmFactor;

  bool mm;

protected:
  DGMesh *mesh;
  INSData *data;

private:
  void cycle(int p);

  double maxEigenValue();
  void setRandomVector(op_dat vec);

  void copy_vec_to_dat(op_dat dat, const double *dat_d);
  void copy_dat_to_vec(op_dat dat, double *dat_d);
  void create_shell_mat(Mat *mat);
  int get_local_unknowns();
  void create_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat);
  void store_vec(Vec *v, op_dat v_dat);

  op_dat bc_dat;

  op_dat tmp_dat[DG_ORDER], u_dat[DG_ORDER], b_dat[DG_ORDER], u_rhs, rhs_rhs;
  double *tmp_dat_data[DG_ORDER], *u_dat_data[DG_ORDER], *b_dat_data[DG_ORDER];
  double *u_rhs_data, *rhs_rhs_data;
  double *factor_data, *mmFactor_data;

  PoissonMat *pMatrix;
};

class PressurePMultigrid : public PMultigrid {
public:
  PressurePMultigrid(DGMesh *m, INSData *d);

  void setup();
};

#endif
