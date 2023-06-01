#ifndef __INS_POISSON_COARSE_MATRIX_H
#define __INS_POISSON_COARSE_MATRIX_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "petscvec.h"
#include "petscksp.h"

#include "dg_mesh/dg_mesh.h"
#include "poisson_matrix.h"

#ifdef INS_BUILD_WITH_AMGX
#include <amgx_c.h>
#endif
#ifdef INS_BUILD_WITH_HYPRE
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#endif

class PoissonCoarseMatrix : public PoissonMatrix {
public:
  virtual void mult(op_dat in, op_dat out) override;
  virtual void multJacobi(op_dat in, op_dat out) override;
  virtual int getUnknowns() override;

  #ifdef INS_BUILD_WITH_AMGX
  bool getAmgXMat(AMGX_matrix_handle** mat);
  #endif
  #ifdef INS_BUILD_WITH_HYPRE
  bool getHYPREMat(HYPRE_ParCSRMatrix** mat);
  void getHYPRERanges(int *ilower, int *iupper, int *jlower, int *jupper);
  #endif

protected:
  virtual void set_glb_ind() override;
  virtual void setPETScMatrix() override;

  #ifdef INS_BUILD_WITH_AMGX
  virtual void setAmgXMatrix();

  AMGX_matrix_handle amgx_mat;
  bool amgx_mat_init = false;
  #endif

  #ifdef INS_BUILD_WITH_HYPRE
  virtual void setHYPREMatrix();

  HYPRE_IJMatrix hypre_mat;
  bool hypre_mat_init = false;
  HYPRE_ParCSRMatrix hypre_parcsr_mat;
  #endif
};

#endif
