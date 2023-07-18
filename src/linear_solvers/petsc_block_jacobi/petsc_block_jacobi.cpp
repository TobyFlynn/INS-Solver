#include "linear_solvers/petsc_block_jacobi.h"

#include "op_seq.h"

#include "utils.h"
#include "linear_solvers/petsc_utils.h"
#include "timing.h"
#include "config.h"
#include "dg_dat_pool.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <type_traits>
#include <stdexcept>

extern Timing *timer;
extern Config *config;
extern DGDatPool *dg_dat_pool;

PETScBlockJacobiSolver::PETScBlockJacobiSolver(DGMesh *m) {
  bc = nullptr;
  nullspace = false;
  pMatInit = false;
  mesh = m;

  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, (DG_FP *)NULL, "block_jacobi_pre");

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  double r_tol, a_tol;
  if(std::is_same<DG_FP,double>::value) {
    r_tol = 1e-8;
    a_tol = 1e-9;
  } else {
    r_tol = 1e-5;
    a_tol = 1e-6;
  }
  config->getDouble("top-level-linear-solvers", "r_tol", r_tol);
  config->getDouble("top-level-linear-solvers", "a_tol", a_tol);
  KSPSetTolerances(ksp, r_tol, a_tol, 1e5, 5e2);
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
  DGTempDat tmp_in  = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_out = dg_dat_pool->requestTempDatCells(DG_NP);
  PETScUtils::copy_vec_to_dat_p_adapt(tmp_in.dat, in_d, mesh);

  matrix->mult(tmp_in.dat, tmp_out.dat);

  PETScUtils::copy_dat_to_vec_p_adapt(tmp_out.dat, out_d, mesh);
  dg_dat_pool->releaseTempDatCells(tmp_in);
  dg_dat_pool->releaseTempDatCells(tmp_out);
  timer->endTimer("PETScBlockJacobiSolver - calc_rhs");
}

// Matrix-free block-jacobi preconditioning function
void PETScBlockJacobiSolver::precond(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScBlockJacobiSolver - precond");
  DGTempDat tmp_in  = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_out = dg_dat_pool->requestTempDatCells(DG_NP);
  PETScUtils::copy_vec_to_dat_p_adapt(tmp_in.dat, in_d, mesh);

  op_par_loop(block_jacobi_pre, "block_jacobi_pre", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(tmp_in.dat,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pre, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_out.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  PETScUtils::copy_dat_to_vec_p_adapt(tmp_out.dat, out_d, mesh);
  dg_dat_pool->releaseTempDatCells(tmp_in);
  dg_dat_pool->releaseTempDatCells(tmp_out);
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
