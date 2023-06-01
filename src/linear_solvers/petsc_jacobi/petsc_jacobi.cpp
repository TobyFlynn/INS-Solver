#include "linear_solvers/petsc_jacobi.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

#include "utils.h"
#include "linear_solvers/petsc_utils.h"
#include "timing.h"
#include "config.h"
#include "dg_dat_pool.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <type_traits>

extern DGConstants *constants;
extern Timing *timer;
extern DGDatPool *dg_dat_pool;
extern Config *config;

void custom_kernel_petsc_pre_jacobi(const int order, char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

PETScJacobiSolver::PETScJacobiSolver(DGMesh *m) {
  bc = nullptr;
  nullspace = false;
  pMatInit = false;
  mesh = m;

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  // KSPSetType(ksp, KSPGMRES);
  KSPSetType(ksp, KSPCG);
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

PETScJacobiSolver::~PETScJacobiSolver() {
  if(pMatInit)
    MatDestroy(&pMat);
  KSPDestroy(&ksp);
  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);
}

void PETScJacobiSolver::init() {
  PETScUtils::create_vec(&b, mesh->cells);
  PETScUtils::create_vec(&x, mesh->cells);
  create_shell_mat();
}

bool PETScJacobiSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("PETScJacobiSolver - solve");

  if(dynamic_cast<PoissonMatrixFreeDiag*>(matrix) == nullptr) {
    throw std::runtime_error("PETScJacobiSolver matrix should be of type PoissonMatrixFreeDiag\n");
  }
  diagMat = dynamic_cast<PoissonMatrixFreeDiag*>(matrix);

  if(nullspace) {
    MatNullSpace ns;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &ns);
    MatSetNullSpace(pMat, ns);
    MatSetTransposeNullSpace(pMat, ns);
    MatNullSpaceDestroy(&ns);
  }

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
  KSPSetOperators(ksp, pMat, pMat);

  if(bc)
    matrix->apply_bc(rhs, bc);

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

  timer->endTimer("PETScJacobiSolver - solve");

  return converged;
}

void PETScJacobiSolver::calc_rhs(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScJacobiSolver - calc_rhs");
  // Copy u to OP2 dat
  DGTempDat tmp_in  = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_out = dg_dat_pool->requestTempDatCells(DG_NP);
  // PETScUtils::copy_vec_to_dat_p_adapt(tmp_in.dat, in_d, mesh);
  PETScUtils::copy_vec_to_dat(tmp_in.dat, in_d);

  matrix->mult(tmp_in.dat, tmp_out.dat);

  // PETScUtils::copy_dat_to_vec_p_adapt(tmp_out.dat, out_d, mesh);
  PETScUtils::copy_dat_to_vec(tmp_out.dat, out_d);
  dg_dat_pool->releaseTempDatCells(tmp_in);
  dg_dat_pool->releaseTempDatCells(tmp_out);
  timer->endTimer("PETScJacobiSolver - calc_rhs");
}

// Matrix-free inv Mass preconditioning function
void PETScJacobiSolver::precond(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScJacobiSolver - precond");
  DGTempDat tmp_in  = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_out = dg_dat_pool->requestTempDatCells(DG_NP);
  // PETScUtils::copy_vec_to_dat_p_adapt(tmp_in.dat, in_d, mesh);
  PETScUtils::copy_vec_to_dat(tmp_in.dat, in_d);

  #if defined(OP2_DG_CUDA) && !defined(USE_OP2_KERNELS)
  custom_kernel_petsc_pre_jacobi(DG_ORDER, "petsc_pre_jacobi", mesh->cells,
              op_arg_dat(diagMat->diag, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_in.dat,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_out.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  #else
  op_par_loop(petsc_pre_jacobi, "petsc_pre_jacobi", mesh->cells,
              op_arg_dat(diagMat->diag, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_in.dat,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_out.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  #endif

  // PETScUtils::copy_dat_to_vec_p_adapt(tmp_out.dat, out_d, mesh);
  PETScUtils::copy_dat_to_vec(tmp_out.dat, out_d);
  dg_dat_pool->releaseTempDatCells(tmp_in);
  dg_dat_pool->releaseTempDatCells(tmp_out);
  timer->endTimer("PETScJacobiSolver - precond");
}
