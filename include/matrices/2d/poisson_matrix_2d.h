#ifndef __INS_POISSON_MATRIX_2D_H
#define __INS_POISSON_MATRIX_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "../poisson_matrix.h"

class PoissonMatrix2D : public PoissonMatrix {
public:
  PoissonMatrix2D(DGMesh2D *m);
  ~PoissonMatrix2D();

  virtual void calc_mat(op_dat bc_types) override;
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;
  virtual void multJacobi(op_dat in, op_dat out) override;
  bool getPETScMat(Mat** mat) override;

  // OP2 Dats
  op_dat h;
  op_dat glb_indBC;
  op_dat orderL, orderR, orderBC;

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc(op_dat bc_types) override;
  virtual void set_glb_ind() override;
  virtual void calc_glb_ind() override;

  DGMesh2D *mesh;

private:
  void setPETScMatrix();
};

#endif
