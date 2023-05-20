#ifndef __AMGX_AMG_H
#define __AMGX_AMG_H

#include "op_seq.h"
#include "linear_solver.h"
#include "dg_mesh/dg_mesh.h"

#ifdef INS_MPI
#include "mpi.h"
#endif

#include <amgx_c.h>

class AmgXAMGSolver : public LinearSolver {
public:
  AmgXAMGSolver(DGMesh *m);
  ~AmgXAMGSolver();

  virtual bool solve(op_dat rhs, op_dat ans) override;

protected:
  DGMesh *mesh;
  KSP ksp;

  bool pMatInit;

  AMGX_matrix_handle *amgx_mat;
  AMGX_vector_handle rhs_amgx;
  AMGX_vector_handle soln_amgx;
  AMGX_solver_handle solver_amgx;

  #ifdef INS_MPI
  MPI_Comm amgx_comm;
  #endif
};

#endif
