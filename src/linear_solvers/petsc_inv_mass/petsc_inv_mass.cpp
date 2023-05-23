#include "linear_solvers/petsc_inv_mass.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

#include "utils.h"
#include "linear_solvers/petsc_utils.h"
#include "timing.h"
#include "dg_dat_pool.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <type_traits>

extern DGConstants *constants;
extern Timing *timer;
extern DGDatPool3D *dg_dat_pool;

PETScInvMassSolver::PETScInvMassSolver(DGMesh *m) {
  bc = nullptr;
  nullspace = false;
  pMatInit = false;
  mesh = m;
  factor = 1.0;
  dat_factor = false;

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
  // PETScUtils::create_vec_p_adapt(&b, matrix->getUnknowns());
  // PETScUtils::create_vec_p_adapt(&x, matrix->getUnknowns());
  PETScUtils::create_vec(&b, mesh->cells);
  PETScUtils::create_vec(&x, mesh->cells);

  // PETScUtils::load_vec_p_adapt(&b, rhs, mesh);
  // PETScUtils::load_vec_p_adapt(&x, ans, mesh);
  PETScUtils::load_vec(&b, rhs);
  PETScUtils::load_vec(&x, ans);

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
  // PETScUtils::store_vec_p_adapt(&solution, ans, mesh);
  PETScUtils::store_vec(&solution, ans);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  timer->endTimer("PETScInvMassSolver - solve");

  return converged;
}

void PETScInvMassSolver::calc_rhs(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScInvMassSolver - calc_rhs");
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
  timer->endTimer("PETScInvMassSolver - calc_rhs");
}

// Matrix-free inv Mass preconditioning function
void PETScInvMassSolver::precond(const DG_FP *in_d, DG_FP *out_d) {
  timer->startTimer("PETScInvMassSolver - precond");
  DGTempDat tmp_in  = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_out = dg_dat_pool->requestTempDatCells(DG_NP);
  // PETScUtils::copy_vec_to_dat_p_adapt(tmp_in.dat, in_d, mesh);
  PETScUtils::copy_vec_to_dat(tmp_in.dat, in_d);

  #if DG_DIM == 3
  if(dat_factor) {
    op_par_loop(petsc_pre_inv_mass_dat, "petsc_pre_inv_mass_dat", mesh->cells,
                op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(factor_dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_in.dat,      -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_out.dat,     -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  } else {
    op_par_loop(petsc_pre_inv_mass, "petsc_pre_inv_mass", mesh->cells,
                op_arg_gbl(&factor,  1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_in.dat,      -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_out.dat,     -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }
  #else
  op_par_loop(petsc_pre_inv_mass_2d, "petsc_pre_inv_mass_2d", mesh->cells,
              op_arg_gbl(constants->get_mat_ptr(DGConstants::INV_MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(&factor,  1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_in.dat,      -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_out.dat,     -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  #endif

  // PETScUtils::copy_dat_to_vec_p_adapt(tmp_out.dat, out_d, mesh);
  PETScUtils::copy_dat_to_vec(tmp_out.dat, out_d);
  dg_dat_pool->releaseTempDatCells(tmp_in);
  dg_dat_pool->releaseTempDatCells(tmp_out);
  timer->endTimer("PETScInvMassSolver - precond");
}

void PETScInvMassSolver::setFactor(const double f) {
  factor = f;
}

void PETScInvMassSolver::setFactor(op_dat f) {
  factor_dat = f;
  dat_factor = true;
}
