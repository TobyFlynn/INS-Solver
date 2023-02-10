#include "linear_solvers/petsc_pmultigrid.h"

#include "op_seq.h"

#include <iostream>

#include "linear_solvers/petsc_utils.h"

PETScPMultigrid::PETScPMultigrid(DGMesh *m) {
  mesh = m;
  nullspace = false;
  pMatInit = false;

  double *tmp_np = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  in  = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "block_jacobi_in");
  out = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "block_jacobi_out");
  free(tmp_np);

  pmultigridSolver = new PMultigridPoissonSolver(mesh);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 2.5e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCSHELL);
  set_shell_pc(pc);
}

PETScPMultigrid::~PETScPMultigrid() {
  MatDestroy(&pMat);
  KSPDestroy(&ksp);
  delete pmultigridSolver;
}

bool PETScPMultigrid::solve(op_dat rhs, op_dat ans) {
  create_shell_mat();
  KSPSetOperators(ksp, pMat, pMat);
  if(nullspace) {
    MatNullSpace ns;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &ns);
    MatSetNullSpace(pMat, ns);
    MatSetTransposeNullSpace(pMat, ns);
    MatNullSpaceDestroy(&ns);
  }

  pmultigridSolver->set_matrix(matrix);
  pmultigridSolver->set_bcs(bc);
  pmultigridSolver->set_nullspace(nullspace);

  if(bc)
    matrix->apply_bc(rhs, bc);

  Vec b, x;
  PETScUtils::create_vec_p_adapt(&b, matrix->unknowns);
  PETScUtils::create_vec_p_adapt(&x, matrix->unknowns);

  PETScUtils::load_vec_p_adapt(&b, rhs, mesh);
  PETScUtils::load_vec_p_adapt(&x, ans, mesh);

  KSPSolve(ksp, b, x);

  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  // Check that the solver converged
  bool converged = true;
  if(reason < 0) {
    double residual;
    KSPGetResidualNorm(ksp, &residual);
    converged = false;
    std::cout << "Number of iterations for linear solver: " << numIt << std::endl;
    std::cout << "Converged reason: " << reason << " Residual: " << residual << std::endl;
  }

  Vec solution;
  KSPGetSolution(ksp, &solution);
  PETScUtils::store_vec_p_adapt(&solution, ans, mesh);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  return converged;
}

void PETScPMultigrid::calc_rhs(const double *in_d, double *out_d) {
  // Copy u to OP2 dat
  PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);

  matrix->mult(in, out);

  PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
}

void PETScPMultigrid::precond(const double *in_d, double *out_d) {
  PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);

  pmultigridSolver->solve(in, out);

  PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
}
