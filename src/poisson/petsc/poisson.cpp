#include "petsc_poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "dg_blas_calls.h"

using namespace std;

extern Timing *timer;

PetscPoissonSolve::PetscPoissonSolve(DGMesh *m, INSData *nsData, LS *s) {
  mesh = m;
  data = nsData;
  ls = s;

  numberIter = 0;
  solveCount = 0;
  pMatInit = false;
  block_jacobi_pre = false;

  u_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rhs_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  in_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  out_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  pre_data = (double *)calloc(DG_NP * DG_NP * mesh->numCells, sizeof(double));

  factor_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  gFactor_data  = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  cFactor_data  = (double *)calloc(DG_CUB_NP * mesh->numCells, sizeof(double));
  mmFactor_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  h_data        = (double *)calloc(mesh->numCells, sizeof(double));
  gDelta_data   = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));

  u   = op_decl_dat(mesh->cells, DG_NP, "double", u_data, "poisson_u");
  rhs = op_decl_dat(mesh->cells, DG_NP, "double", rhs_data, "poisson_rhs");
  in  = op_decl_dat(mesh->cells, DG_NP, "double", in_data, "poisson_in");
  out = op_decl_dat(mesh->cells, DG_NP, "double", out_data, "poisson_out");
  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", pre_data, "poisson_pre");

  factor   = op_decl_dat(mesh->cells, DG_NP, "double", factor_data, "poisson_factor");
  gFactor  = op_decl_dat(mesh->cells, DG_G_NP, "double", gFactor_data, "poisson_gFactor");
  cFactor  = op_decl_dat(mesh->cells, DG_CUB_NP, "double", cFactor_data, "poisson_cFactor");
  mmFactor = op_decl_dat(mesh->cells, DG_NP, "double", mmFactor_data, "poisson_mmFactor");
  h        = op_decl_dat(mesh->cells, 1, "double", h_data, "poisson_h");
  gDelta   = op_decl_dat(mesh->cells, DG_G_NP, "double", gDelta_data, "poisson_gDelta");

  mat = new PoissonMat(mesh);
}

PetscPoissonSolve::~PetscPoissonSolve() {
  free(u_data);
  free(rhs_data);
  free(in_data);
  free(out_data);
  free(pre_data);

  free(factor_data);
  free(gFactor_data);
  free(cFactor_data);
  free(mmFactor_data);
  free(h_data);
  free(gDelta_data);

  if(pMatInit)
    MatDestroy(&pMat);
  
  delete mat;
}

void PetscPoissonSolve::init() {
  mat->init();
}

bool PetscPoissonSolve::solve(op_dat b_dat, op_dat x_dat) {
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 1e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

  KSPSetOperators(ksp, pMat, pMat);

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

  op_par_loop(poisson_apply_bc, "poisson_apply_bc", mesh->bedges,
              op_arg_dat(mesh->order,     0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mat->op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(b_dat,  0, mesh->bedge2cells, DG_NP, "double", OP_INC));

  create_vec(&b);
  create_vec(&x);

  load_vec(&b, b_dat);
  load_vec(&x, x_dat);

  timer->startTimer("KSPSolve");
  KSPSolve(ksp, b, x);
  timer->endTimer("KSPSolve");

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

  destroy_vec(&b);
  destroy_vec(&x);

  KSPDestroy(&ksp);

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

PetscPressureSolve::PetscPressureSolve(DGMesh *m, INSData *d, LS *s) : PetscPoissonSolve(m, d, s) {}

PetscViscositySolve::PetscViscositySolve(DGMesh *m, INSData *d, LS *s) : PetscPoissonSolve(m, d, s) {}

void PetscPressureSolve::setup() {
  mat->update_glb_ind();

  massMat = false;

  op_par_loop(poisson_pr_fact, "poisson_pr_fact", mesh->cells,
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE));

  timer->startTimer("Build Pr Mat");
  mat->calc_mat(factor);
  setMatrix();
  timer->endTimer("Build Pr Mat");
}

void PetscViscositySolve::setup(double mmConst) {
  mat->update_glb_ind();

  massMat = true;
  massFactor = mmConst;
  block_jacobi_pre = true;

  op_par_loop(poisson_vis_fact, "poisson_vis_fact", mesh->cells,
              op_arg_gbl(&mmConst,   1, "double", OP_READ),
              op_arg_dat(data->mu,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(mmFactor,  -1, OP_ID, DG_NP, "double", OP_WRITE));

  timer->startTimer("Build Vis Mat");
  mat->calc_mat_mm(factor, mmFactor);

  if(block_jacobi_pre) {
    inv_blas(mesh, mat->op1, pre);
  }

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
