#ifndef __INS_CUB_POISSON_MATRIX_2D_H
#define __INS_CUB_POISSON_MATRIX_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "poisson_matrix_over_int_2d.h"

class CubPoissonMatrix2D : public PoissonMatrixOverInt2D {
public:
  CubPoissonMatrix2D(DGMesh2D *m);

protected:
  virtual void calc_op1() override;
};

#endif
