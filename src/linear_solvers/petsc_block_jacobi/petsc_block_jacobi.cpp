#include "linear_solvers/petsc_block_jacobi.h"

#include "op_seq.h"

#include "utils.h"
#include "linear_solvers/petsc_utils.h"
#include "timing.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <type_traits>

extern Timing *timer;

PETScBlockJacobiSolver::PETScBlockJacobiSolver(DGMesh *m) {
  bc = nullptr;
  nullspace = false;
  pMatInit = false;
  mesh = m;

  DG_FP *tmp_np = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *tmp_np_np = (DG_FP *)calloc(DG_NP * DG_NP * mesh->cells->size, sizeof(DG_FP));

  in  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "block_jacobi_in");
  out = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "block_jacobi_out");
  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, tmp_np_np, "block_jacobi_pre");

  free(tmp_np_np);
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

PETScBlockJacobiSolver::~PETScBlockJacobiSolver() {
  MatDestroy(&pMat);
  KSPDestroy(&ksp);
}

bool PETScBlockJacobiSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("PETScBlockJacobiSolver - solve");
  // TODO only call when necessary
  calc_precond_mat();
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

  timer->startTimer("PETScBlockJacobiSolver - KSPSolve");
  KSPSolve(ksp, b, x);
  timer->endTimer("PETScBlockJacobiSolver - KSPSolve");

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

  timer->endTimer("PETScBlockJacobiSolver - solve");

  return converged;
}

void PETScBlockJacobiSolver::calc_rhs(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScBlockJacobiSolver - calc_rhs");
  // Copy u to OP2 dat
  PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);

  matrix->mult(in, out);

  PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
  timer->endTimer("PETScBlockJacobiSolver - calc_rhs");
}

// Matrix-free block-jacobi preconditioning function
void PETScBlockJacobiSolver::precond(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScBlockJacobiSolver - precond");
  PETScUtils::copy_vec_to_dat_p_adapt(in, in_d, mesh);

  op_par_loop(block_jacobi_pre, "block_jacobi_pre", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pre, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  PETScUtils::copy_dat_to_vec_p_adapt(out, out_d, mesh);
  timer->endTimer("PETScBlockJacobiSolver - precond");
}

void PETScBlockJacobiSolver::calc_precond_mat() {
  timer->startTimer("PETScBlockJacobiSolver - calc_precond_mat");
  const DG_FP *op1_ptr = getOP2PtrHost(matrix->op1, OP_READ);
  DG_FP *pre_ptr = getOP2PtrHost(pre, OP_WRITE);

  #pragma omp parallel for
  for(int i = 0; i < mesh->cells->size; i++) {
    const DG_FP *in_c = op1_ptr + i * matrix->op1->dim;
    DG_FP *inv_c      = pre_ptr + i * pre->dim;

    arma::Mat<DG_FP> a(in_c, DG_NP, DG_NP);
    arma::Mat<DG_FP> b(inv_c, DG_NP, DG_NP, false, true);

    #ifdef DG_COL_MAJ
    b = arma::inv(a);
    #else
    b = arma::inv(a.t()).t();
    #endif
    // b = arma::inv_sympd(a);
  }

  releaseOP2PtrHost(matrix->op1, OP_READ, op1_ptr);
  releaseOP2PtrHost(pre, OP_WRITE, pre_ptr);
  timer->endTimer("PETScBlockJacobiSolver - calc_precond_mat");
}
