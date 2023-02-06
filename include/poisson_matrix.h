#ifndef __INS_POISSON_MATRIX_H
#define __INS_POISSON_MATRIX_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"

class PoissonMatrix2D {
public:
  PoissonMatrix2D(DGMesh2D *m);

  void init();
  void calc_mat(op_dat fact);
  void calc_mat_mm(op_dat fact, op_dat mmFact);
  void mult(op_dat in, op_dat out);
  void multJacobi(op_dat in, op_dat out);
  void transpose();

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);

  void update_glb_ind();

  // OP2 Dats
  op_dat op1, op2[2], op_bc, h;
  op_dat factor, mmFactor, cFactor, gFactor;
  op_dat glb_ind, glb_indL, glb_indR, glb_indBC;
  op_dat orderL, orderR, orderBC;

  int unknowns;

private:
  void calc_cub_sub_mat();
  void calc_gauss_sub_mat();
  void calc_mm_mat();

  int get_local_unknowns();
  void set_glb_ind();

  DGMesh2D *mesh;

  int dirichlet[3];
  int neumann[3];
};

#endif
