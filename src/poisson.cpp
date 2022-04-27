#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "dg_blas_calls.h"
#include "petsc_utils.h"

using namespace std;

extern Timing *timer;

PoissonSolve::PoissonSolve(DGMesh *m, INSData *nsData, LS *s) {
  mesh = m;
  data = nsData;
  ls = s;

  pMatrix = new PoissonMat(mesh);

  numberIter = 0;
  solveCount = 0;
  pMatInit = false;
  block_jacobi_pre = false;

  glb_ind_data   = (int *)calloc(mesh->numCells, sizeof(int));
  glb_indL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indBC_data = (int *)calloc(mesh->numBoundaryEdges, sizeof(int));

  orderL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  orderR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  orderBC_data = (int *)calloc(mesh->numEdges, sizeof(int));

  u_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rhs_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  in_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  out_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  pre_data = (double *)calloc(DG_NP * DG_NP * mesh->numCells, sizeof(double));

  factor_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  mmFactor_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  glb_ind   = op_decl_dat(mesh->cells, 1, "int", glb_ind_data, "poisson_glb_ind");
  glb_indL  = op_decl_dat(mesh->edges, 1, "int", glb_indL_data, "poisson_glb_indL");
  glb_indR  = op_decl_dat(mesh->edges, 1, "int", glb_indR_data, "poisson_glb_indR");
  glb_indBC = op_decl_dat(mesh->bedges, 1, "int", glb_indBC_data, "poisson_glb_indBC");

  orderL  = op_decl_dat(mesh->edges, 1, "int", orderL_data, "poisson_orderL");
  orderR  = op_decl_dat(mesh->edges, 1, "int", orderR_data, "poisson_orderR");
  orderBC = op_decl_dat(mesh->bedges, 1, "int", orderBC_data, "poisson_orderBC");

  u   = op_decl_dat(mesh->cells, DG_NP, "double", u_data, "poisson_u");
  rhs = op_decl_dat(mesh->cells, DG_NP, "double", rhs_data, "poisson_rhs");
  in  = op_decl_dat(mesh->cells, DG_NP, "double", in_data, "poisson_in");
  out = op_decl_dat(mesh->cells, DG_NP, "double", out_data, "poisson_out");
  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", pre_data, "poisson_pre");

  factor   = op_decl_dat(mesh->cells, DG_NP, "double", factor_data, "poisson_factor");
  mmFactor = op_decl_dat(mesh->cells, DG_NP, "double", mmFactor_data, "poisson_mmFactor");
}

PoissonSolve::~PoissonSolve() {
  free(glb_ind_data);
  free(glb_indL_data);
  free(glb_indR_data);
  free(glb_indBC_data);

  free(orderL_data);
  free(orderR_data);
  free(orderBC_data);

  free(u_data);
  free(rhs_data);
  free(in_data);
  free(out_data);
  free(pre_data);

  free(factor_data);
  free(mmFactor_data);

  delete pMatrix;

  if(pMatInit)
    MatDestroy(&pMat);
}

void PoissonSolve::init() {
  pMatrix->init();
}

void PoissonSolve::update_glb_ind() {
  setGlbInd();
  op_par_loop(copy_to_edges, "copy_to_edges", mesh->edges,
              op_arg_dat(glb_ind, -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));
  op_par_loop(copy_to_bedges, "copy_to_bedges", mesh->bedges,
              op_arg_dat(glb_ind, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indBC, -1, OP_ID, 1, "int", OP_WRITE));

  op_par_loop(copy_to_edges, "copy_to_edges", mesh->edges,
              op_arg_dat(mesh->order, -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_WRITE));
  op_par_loop(copy_to_bedges, "copy_to_bedges", mesh->bedges,
              op_arg_dat(mesh->order, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(orderBC, -1, OP_ID, 1, "int", OP_WRITE));
}

bool PoissonSolve::solve(op_dat b_dat, op_dat x_dat) {
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
              op_arg_dat(pMatrix->op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(b_dat,  0, mesh->bedge2cells, DG_NP, "double", OP_INC));

  PETScUtils::create_vec(&b, unknowns);
  PETScUtils::create_vec(&x, unknowns);

  PETScUtils::dat_to_vec(&b, b_dat);
  PETScUtils::dat_to_vec(&x, x_dat);

  timer->startKSPSolve();
  KSPSolve(ksp, b, x);
  timer->endKSPSolve();

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
  PETScUtils::vec_to_dat(&solution, x_dat);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  KSPDestroy(&ksp);

  return converged;
}

// Matrix-free Mat-Vec mult function
void PoissonSolve::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  PETScUtils::ptr_to_dat(u, u_d);

  pMatrix->mult(u, rhs);

  PETScUtils::dat_to_ptr(rhs, rhs_d);
}

// Matrix-free block-jacobi preconditioning function
void PoissonSolve::precond(const double *in_d, double *out_d) {
  PETScUtils::ptr_to_dat(in, in_d);

  op_par_loop(poisson_pre, "poisson_pre", mesh->cells,
              op_arg_dat(in,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(pre, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_WRITE));

  PETScUtils::dat_to_ptr(out, out_d);
}

void PoissonSolve::setDirichletBCs(int *d) {
  pMatrix->setDirichletBCs(d);
}

void PoissonSolve::setNeumannBCs(int *n) {
  pMatrix->setNeumannBCs(n);
}

void PoissonSolve::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double PoissonSolve::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}

PressureSolve::PressureSolve(DGMesh *m, INSData *d, LS *s) : PoissonSolve(m, d, s) {}

ViscositySolve::ViscositySolve(DGMesh *m, INSData *d, LS *s) : PoissonSolve(m, d, s) {}

void PressureSolve::setup() {
  unknowns = mesh->get_local_vec_unknowns();
  update_glb_ind();

  massMat = false;

  op_par_loop(poisson_pr_fact, "poisson_pr_fact", mesh->cells,
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE));

  pMatrix->calc_mat(factor);

  timer->startBuildMat();
  setMatrix();
  timer->endBuildMat();
}

void ViscositySolve::setup(double mmConst) {
  unknowns = mesh->get_local_vec_unknowns();
  update_glb_ind();

  massMat = true;
  massFactor = mmConst;
  block_jacobi_pre = true;

  op_par_loop(poisson_vis_fact, "poisson_vis_fact", mesh->cells,
              op_arg_gbl(&mmConst,   1, "double", OP_READ),
              op_arg_dat(data->mu,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(mmFactor,  -1, OP_ID, DG_NP, "double", OP_WRITE));

  pMatrix->calc_mat_mm(factor, mmFactor);

  if(block_jacobi_pre) {
    inv_blas(mesh, pMatrix->op1, pre);
  }

  create_shell_mat();
}
