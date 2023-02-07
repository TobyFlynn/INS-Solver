#include "linear_solvers/petsc_poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "timing.h"
#include "utils.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

using namespace std;

extern Timing *timer;

PetscPoissonSolve::PetscPoissonSolve(DGMesh2D *m) {
  mesh = m;

  numberIter = 0;
  solveCount = 0;
  pMatInit = false;
  block_jacobi_pre = false;
  vec_created = false;
  massMat = false;
  prev_unknowns = 0;

  double *tmp_np = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  double *tmp_np_np = (double *)calloc(DG_NP * DG_NP * mesh->cells->size, sizeof(double));
  double *tmp_g_np = (double *)calloc(DG_G_NP * mesh->cells->size, sizeof(double));
  double *tmp_cub_np = (double *)calloc(DG_CUB_NP * mesh->cells->size, sizeof(double));
  double *tmp_1 = (double *)calloc(mesh->cells->size, sizeof(double));

  u   = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "poisson_u");
  rhs = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "poisson_rhs");
  in  = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "poisson_in");
  out = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "poisson_out");
  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", tmp_np_np, "poisson_pre");

  factor   = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "poisson_factor");
  gFactor  = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_g_np, "poisson_gFactor");
  cFactor  = op_decl_dat(mesh->cells, DG_CUB_NP, "double", tmp_cub_np, "poisson_cFactor");
  mmFactor = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "poisson_mmFactor");
  h        = op_decl_dat(mesh->cells, 1, "double", tmp_1, "poisson_h");
  gDelta   = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_g_np, "poisson_gDelta");

  free(tmp_1);
  free(tmp_cub_np);
  free(tmp_g_np);
  free(tmp_np_np);
  free(tmp_np);

  mat = new PoissonMatrix2D(mesh);
}

PetscPoissonSolve::~PetscPoissonSolve() {
  if(pMatInit) {
    KSPDestroy(&ksp);
    MatDestroy(&pMat);
  }
  if(vec_created) {
    destroy_vec(&b);
    destroy_vec(&x);
  }

  delete mat;
}

void PetscPoissonSolve::init() {
  mat->init();

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 1e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

  PC pc;
  KSPGetPC(ksp, &pc);

  if(massMat) {
    PCSetType(pc, PCSHELL);
    set_shell_pc(pc);
  } else {
    PCSetType(pc, PCGAMG);
    PCGAMGSetNSmooths(pc, 4);
    PCGAMGSetSquareGraph(pc, 1);
    PCGAMGSetNlevels(pc, 20);
    PCMGSetLevels(pc, 20, NULL);
    PCMGSetCycleType(pc, PC_MG_CYCLE_W);
    PCGAMGSetRepartition(pc, PETSC_TRUE);
    PCGAMGSetReuseInterpolation(pc, PETSC_TRUE);
  }
}

bool PetscPoissonSolve::solve(op_dat b_dat, op_dat x_dat) {
  timer->startTimer("PETSc - solve");
  // if(!vec_created) {
    create_vec(&b);
    create_vec(&x);
  //   vec_created = true;
  // }

  KSPSetOperators(ksp, pMat, pMat);

  if(mesh->bface2cells) {
    op_par_loop(poisson_apply_bc, "poisson_apply_bc", mesh->bfaces,
                op_arg_dat(mesh->order,     0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mat->op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_dat(bc_dat, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(b_dat,  0, mesh->bface2cells, DG_NP, "double", OP_INC));
  }

  load_vec(&b, b_dat);
  load_vec(&x, x_dat);

  timer->startTimer("PETSc - KSPSolve");
  KSPSolve(ksp, b, x);
  timer->endTimer("PETSc - KSPSolve");

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
    cout << "Number of iterations for linear solver: " << numIt << endl;
    cout << "Converged reason: " << reason << " Residual: " << residual << endl;
  }
  numberIter += numIt;
  solveCount++;

  // Get solution
  Vec solution;
  KSPGetSolution(ksp, &solution);
  store_vec(&solution, x_dat);
  timer->endTimer("PETSc - solve");

  destroy_vec(&b);
  destroy_vec(&x);

  return converged;
}

// Matrix-free Mat-Vec mult function
void PetscPoissonSolve::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  copy_vec_to_dat(u, u_d);

  mat->mult(u, rhs);

  copy_dat_to_vec(rhs, rhs_d);
}

// Matrix-free block-jacobi preconditioning function
void PetscPoissonSolve::precond(const double *in_d, double *out_d) {
  copy_vec_to_dat(in, in_d);

  op_par_loop(poisson_pre, "poisson_pre", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(pre, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_WRITE));

  copy_dat_to_vec(out, out_d);
}

void PetscPoissonSolve::setDirichletBCs(int *d) {
  mat->setDirichletBCs(d);
}

void PetscPoissonSolve::setNeumannBCs(int *n) {
  mat->setNeumannBCs(n);
}

void PetscPoissonSolve::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double PetscPoissonSolve::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}

PetscPressureSolve::PetscPressureSolve(DGMesh2D *m) : PetscPoissonSolve(m) {}

PetscViscositySolve::PetscViscositySolve(DGMesh2D *m) : PetscPoissonSolve(m) {
  massMat = true;
}

void PetscPressureSolve::setup(op_dat rho) {
  mat->update_glb_ind();

  massMat = false;

  op_par_loop(poisson_pr_fact, "poisson_pr_fact", mesh->cells,
              op_arg_dat(rho,    -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor, -1, OP_ID, DG_NP, "double", OP_WRITE));

  timer->startTimer("Build Pr Mat");
  mat->calc_mat(factor);
  setMatrix();
  if(!mesh->bface2cells) {
    MatNullSpace nullspace;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
    MatSetNullSpace(pMat, nullspace);
    MatSetTransposeNullSpace(pMat, nullspace);
    MatNullSpaceDestroy(&nullspace);
  }
  timer->endTimer("Build Pr Mat");
}

void PetscViscositySolve::setup(double mmConst, op_dat rho, op_dat mu) {
  mat->update_glb_ind();

  massMat = true;
  massFactor = mmConst;
  block_jacobi_pre = true;

  op_par_loop(poisson_vis_fact, "poisson_vis_fact", mesh->cells,
              op_arg_gbl(&mmConst,   1, "double", OP_READ),
              op_arg_dat(mu,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(mmFactor,  -1, OP_ID, DG_NP, "double", OP_WRITE));

  timer->startTimer("Build Vis Mat");
  mat->calc_mat_mm(factor, mmFactor);

  timer->startTimer("Vis Mat Pre Inv");
  if(block_jacobi_pre) {
    calc_precond_mat();
  }
  timer->endTimer("Vis Mat Pre Inv");

  create_shell_mat();
  timer->endTimer("Build Vis Mat");

  // timer->startBuildMat();
  // setMatrix();
  // timer->endBuildMat();

  // create_shell_mat();

  // PC pc;
  // KSPGetPC(ksp, &pc);
  // PCSetType(pc, PCSHELL);
  // set_shell_pc(pc);
}

void PetscViscositySolve::calc_precond_mat() {
  const double *op1_ptr = getOP2PtrHost(mat->op1, OP_READ);
  double *pre_ptr = getOP2PtrHost(pre, OP_WRITE);

  for(int i = 0; i < mesh->cells->size; i++) {
    const double *in_c = op1_ptr + i * mat->op1->dim;
    double *inv_c      = pre_ptr + i * pre->dim;

    arma::mat a(in_c, DG_NP, DG_NP);
    arma::mat b(inv_c, DG_NP, DG_NP, false, true);

    b = arma::inv(a);
    // b = arma::inv_sympd(a);
  }

  releaseOP2PtrHost(mat->op1, OP_READ, op1_ptr);
  releaseOP2PtrHost(pre, OP_WRITE, pre_ptr);
}
