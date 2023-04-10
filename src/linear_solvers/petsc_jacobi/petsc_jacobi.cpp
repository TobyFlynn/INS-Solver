#include "linear_solvers/petsc_jacobi.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

#include "utils.h"
#include "linear_solvers/petsc_utils.h"
#include "timing.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <type_traits>

extern DGConstants *constants;
extern Timing *timer;

PETScJacobiSolver::PETScJacobiSolver(DGMesh *m) {
  bc = nullptr;
  nullspace = false;
  pMatInit = false;
  mesh = m;

  DG_FP *tmp_np = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));

  in  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "jacobi_in");
  out = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "jacobi_out");

  free(tmp_np);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  if(std::is_same<DG_FP,double>::value)
    KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 5e2);
  else
    KSPSetTolerances(ksp, 1e-5, 1e-50, 1e5, 5e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCSHELL);
  set_shell_pc(pc);
}

PETScJacobiSolver::~PETScJacobiSolver() {
  if(pMatInit)
    MatDestroy(&pMat);
  KSPDestroy(&ksp);
}

bool PETScJacobiSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("PETScJacobiSolver - solve");

  if(dynamic_cast<PoissonMatrixFreeDiag*>(matrix) == nullptr) {
    throw std::runtime_error("PETScJacobiSolver matrix should be of type PoissonMatrixFreeDiag\n");
  }
  diagMat = dynamic_cast<PoissonMatrixFreeDiag*>(matrix);

  create_shell_mat();
  KSPSetOperators(ksp, pMat, pMat);
  if(nullspace) {
    MatNullSpace ns;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &ns);
    MatSetNullSpace(pMat, ns);
    MatSetTransposeNullSpace(pMat, ns);
    MatNullSpaceDestroy(&ns);
  }

  if(bc)
    matrix->apply_bc(rhs, bc);

  Vec b, x;
  PETScUtils::create_vec_p_adapt(&b, matrix->getUnknowns());
  PETScUtils::create_vec_p_adapt(&x, matrix->getUnknowns());

  // PETScUtils::load_vec_p_adapt(&b, rhs, mesh);
  // PETScUtils::load_vec_p_adapt(&x, ans, mesh);
  PETScUtils::load_vec(&b, rhs);
  PETScUtils::load_vec(&x, ans);

  timer->startTimer("PETScJacobiSolver - KSPSolve");
  KSPSolve(ksp, b, x);
  timer->endTimer("PETScJacobiSolver - KSPSolve");

  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  // Check that the solver converged
  bool converged = true;
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

  timer->endTimer("PETScJacobiSolver - solve");

  return converged;
}

void PETScJacobiSolver::calc_rhs(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScJacobiSolver - calc_rhs");
  // Copy u to OP2 dat
  // PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);
  PETScUtils::copy_vec_to_dat(in, in_d);

  matrix->mult(in, out);

  // PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
  PETScUtils::copy_dat_to_vec(out, out_d);
  timer->endTimer("PETScJacobiSolver - calc_rhs");
}

// Matrix-free inv Mass preconditioning function
void PETScJacobiSolver::precond(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScJacobiSolver - precond");
  // PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);
  PETScUtils::copy_vec_to_dat(in, in_d);

  op_par_loop(petsc_pre_jacobi, "petsc_pre_jacobi", mesh->cells,
              op_arg_dat(diagMat->diag, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
  PETScUtils::copy_dat_to_vec(out, out_d);
  timer->endTimer("PETScJacobiSolver - precond");
}
