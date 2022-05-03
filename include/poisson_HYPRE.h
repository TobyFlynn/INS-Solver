#ifndef __POISSON_HYPRE_H
#define __POISSON_HYPRE_H

#include "op_seq.h"
#include "ins_data.h"
#include "timing.h"
#include "dg_mesh.h"
#include "ls.h"
#include "poisson_mat.h"

#include "HYPRE_config.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

class PoissonSolveHYPRE {
public:
  PoissonSolveHYPRE(DGMesh *m, INSData *nsData, LS *s);
  ~PoissonSolveHYPRE();

  void init();
  bool solve(op_dat b_dat, op_dat x_dat);
  /*
  void calc_rhs(const double *u_d, double *rhs_d);
  void precond(const double *in_d, double *out_d);

  double getAverageConvergeIter();
  */

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);

  // OP2 Dats
  op_dat glb_ind, glb_indL, glb_indR, glb_indBC;
  op_dat u, rhs, in, out, pre;
  op_dat factor, mmFactor;
  op_dat orderL, orderR, orderBC;

  int unknowns;

  PoissonMat *pMatrix;

protected:
  void setMatrix();
  /*
  void create_shell_mat();
  void set_shell_pc(PC pc);
  */
  void update_glb_ind();

  DGMesh *mesh;
  INSData *data;
  LS *ls;

  HYPRE_IJMatrix mat;
  HYPRE_ParCSRMatrix parcsr_mat;
  HYPRE_Solver precon, solver;

  HYPRE_IJVector   ij_x, ij_b;
  HYPRE_ParVector  par_x, par_b;

  int dirichlet[3];
  int neumann[3];

private:
  void setGlbInd();

  op_dat bc_dat;

  int numberIter, solveCount;

  int *glb_ind_data, *glb_indL_data, *glb_indR_data, *glb_indBC_data;
  double *u_data, *rhs_data, *in_data, *out_data, *pre_data;
  double *factor_data, *mmFactor_data;
  int *orderL_data, *orderR_data, *orderBC_data;
};

class PressureSolveHYPRE : public PoissonSolveHYPRE {
public:
  PressureSolveHYPRE(DGMesh *m, INSData *d, LS *s);

  void setup();
};

#endif
