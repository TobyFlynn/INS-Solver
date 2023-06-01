#ifndef __HYPRE_AMG_H
#define __HYPRE_AMG_H

#if defined(INS_BUILD_WITH_HYPRE)

#include "op_seq.h"
#include "linear_solver.h"
#include "dg_mesh/dg_mesh.h"

#ifdef INS_MPI
#include "mpi.h"
#endif

#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

class HYPREAMGSolver : public LinearSolver {
public:
  HYPREAMGSolver(DGMesh *m);
  ~HYPREAMGSolver();

  virtual bool solve(op_dat rhs, op_dat ans) override;

protected:
  DGMesh *mesh;
  KSP ksp;

  bool pMatInit;
  bool vec_init;

  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
  HYPRE_Solver solver, precond;

  #ifdef INS_MPI
  MPI_Comm hypre_comm;
  #endif

  int max_pcg_iter, pcg_print_level, pcg_logging, amg_print_level;
  int amg_coarsen_type, amg_relax_type, amg_num_sweeps, amg_iter;
  int amg_keep_transpose, amg_rap2, amg_module_rap2;
  double pcg_tol, amg_strong_threshold, amg_trunc_factor;
};

#endif

#endif
