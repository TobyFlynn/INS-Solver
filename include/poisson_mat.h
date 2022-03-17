#ifndef __INS_POISSON_MAT_H
#define __INS_POISSON_MAT_H

#include "op_seq.h"
#include "dg_mesh.h"

class PoissonMat {
public:
  PoissonMat(DGMesh *m);
  ~PoissonMat();

  void init();
  void calc_mat(op_dat fact);
  void calc_mat_mm(op_dat fact, op_dat mmFact);
  void mult(op_dat in, op_dat out);
  void multJacobi(op_dat in, op_dat out);

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);

  // OP2 Dats
  op_dat op1, op2[2], op_bc, h;
  op_dat factor, mmFactor, cFactor, gFactor;

  int unknowns;

private:
  void calc_cub_sub_mat();
  void calc_gauss_sub_mat();
  void calc_mm_mat();

  DGMesh *mesh;

  int dirichlet[3];
  int neumann[3];

  double *op1_data, *op2_data[2], *op_bc_data, *h_data, *cFactor_data, *gFactor_data;
};

#endif
