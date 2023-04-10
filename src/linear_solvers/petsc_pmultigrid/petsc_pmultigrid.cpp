#include "linear_solvers/petsc_pmultigrid.h"

#include "op_seq.h"

#include <iostream>
#include <type_traits>

#include "linear_solvers/petsc_utils.h"
#include "timing.h"

extern Timing *timer;

PETScPMultigrid::PETScPMultigrid(DGMesh *m) {
  bc = nullptr;
  mesh = m;
  nullspace = false;
  pMatInit = false;

  DG_FP *tmp_np = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  in  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "petsc_pmultigrid_in");
  out = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "petsc_pmultigrid_out");
  free(tmp_np);

  pmultigridSolver = new PMultigridPoissonSolver(mesh);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  if(std::is_same<DG_FP,double>::value)
    KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 2.5e2);
  else
    KSPSetTolerances(ksp, 1e-5, 1e-50, 1e5, 2.5e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCSHELL);
  set_shell_pc(pc);
}

PETScPMultigrid::~PETScPMultigrid() {
  if(pMatInit)
    MatDestroy(&pMat);
  KSPDestroy(&ksp);
  delete pmultigridSolver;
}

void PETScPMultigrid::set_coarse_matrix(PoissonCoarseMatrix *c_mat) {
  pmultigridSolver->set_coarse_matrix(c_mat);
}

bool PETScPMultigrid::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("PETScPMultigrid - solve");
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
  pmultigridSolver->set_nullspace(nullspace);

  if(bc)
    matrix->apply_bc(rhs, bc);

  Vec b, x;
  PETScUtils::create_vec_p_adapt(&b, matrix->getUnknowns());
  PETScUtils::create_vec_p_adapt(&x, matrix->getUnknowns());

  // PETScUtils::load_vec_p_adapt(&b, rhs, mesh);
  // PETScUtils::load_vec_p_adapt(&x, ans, mesh);
  PETScUtils::load_vec(&b, rhs);
  PETScUtils::load_vec(&x, ans);

  timer->startTimer("PETScPMultigrid - KSPSolve");
  KSPSolve(ksp, b, x);
  timer->endTimer("PETScPMultigrid - KSPSolve");

  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  // Check that the solver converged
  bool converged = true;
  // std::cout << "Number of iterations for PETSc PMultigrid: " << numIt << std::endl;
  if(reason < 0) {
    DG_FP residual;
    KSPGetResidualNorm(ksp, &residual);
    converged = false;
    std::cout << "Number of iterations for linear solver: " << numIt << std::endl;
    std::cout << "Converged reason: " << reason << " Residual: " << residual << std::endl;
  }

  Vec solution;
  KSPGetSolution(ksp, &solution);
  // PETScUtils::store_vec_p_adapt(&solution, ans, mesh);
  PETScUtils::store_vec(&solution, ans);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  timer->endTimer("PETScPMultigrid - solve");

  return converged;
}

void PETScPMultigrid::calc_rhs(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScPMultigrid - calc_rhs");
  // Copy u to OP2 dat
  // PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);
  PETScUtils::copy_vec_to_dat(in, in_d);

  matrix->mult(in, out);

  // PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
  PETScUtils::copy_dat_to_vec(out, out_d);
  timer->endTimer("PETScPMultigrid - calc_rhs");
}

void PETScPMultigrid::precond(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScPMultigrid - precond");
  // PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);
  PETScUtils::copy_vec_to_dat(in, in_d);

  pmultigridSolver->solve(in, out);

  // PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
  PETScUtils::copy_dat_to_vec(out, out_d);
  timer->endTimer("PETScPMultigrid - precond");
}
