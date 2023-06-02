#ifndef __INS_POISSON_MATRIX_OVER_INT_2D_H
#define __INS_POISSON_MATRIX_OVER_INT_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "../poisson_matrix.h"

class PoissonMatrixOverInt2D : public PoissonMatrix {
public:
  PoissonMatrixOverInt2D(DGMesh2D *m);
  ~PoissonMatrixOverInt2D();

  virtual void calc_mat() override;
  virtual void apply_bc(op_dat rhs, op_dat bc) override;

  // OP2 Dats
  op_dat h;
  op_dat glb_indBC;
  op_dat orderBC;

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;
  virtual void calc_glb_ind() override;

  DGMesh2D *mesh;
};

#endif
