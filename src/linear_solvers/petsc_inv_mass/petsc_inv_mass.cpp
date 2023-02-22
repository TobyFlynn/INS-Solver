#include "linear_solvers/petsc_inv_mass.h"

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

PETScInvMassSolver::PETScInvMassSolver(DGMesh *m) {
  nullspace = false;
  pMatInit = false;
  mesh = m;
  factor = 1.0;

  DG_FP *tmp_np = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));

  in  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "block_jacobi_in");
  out = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "block_jacobi_out");

  free(tmp_np);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  if(std::is_same<DG_FP,double>::value)
    KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 5e2);
  else
    KSPSetTolerances(ksp, 1e-6, 1e-50, 1e5, 5e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCSHELL);
  set_shell_pc(pc);
}

PETScInvMassSolver::~PETScInvMassSolver() {
  MatDestroy(&pMat);
  KSPDestroy(&ksp);
}

bool PETScInvMassSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("PETScInvMassSolver - solve");
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

  PETScUtils::load_vec_p_adapt(&b, rhs, mesh);
  PETScUtils::load_vec_p_adapt(&x, ans, mesh);

  timer->startTimer("PETScInvMassSolver - KSPSolve");
  KSPSolve(ksp, b, x);
  timer->endTimer("PETScInvMassSolver - KSPSolve");

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
  PETScUtils::store_vec_p_adapt(&solution, ans, mesh);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  timer->endTimer("PETScInvMassSolver - solve");

  return converged;
}

void PETScInvMassSolver::calc_rhs(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScInvMassSolver - calc_rhs");
  // Copy u to OP2 dat
  PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);

  matrix->mult(in, out);

  PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
  timer->endTimer("PETScInvMassSolver - calc_rhs");
}

// Matrix-free inv Mass preconditioning function
void PETScInvMassSolver::precond(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScInvMassSolver - precond");
  PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);

  op_par_loop(petsc_pre_inv_mass, "petsc_pre_inv_mass", mesh->cells,
              op_arg_gbl(constants->get_mat_ptr(DGConstants::INV_MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(&factor,  1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(in,      -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out,     -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
  timer->endTimer("PETScInvMassSolver - precond");
}

void PETScInvMassSolver::setFactor(const double f) {
  factor = f;
}
