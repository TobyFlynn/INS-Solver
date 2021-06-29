#include "viscosity_solve.h"

#include "op_seq.h"

#include <iostream>

#include "constants.h"
#include "blas_calls.h"

extern Constants *constants;

using namespace std;

ViscositySolve::ViscositySolve(INSData *d, CubatureData *c, GaussData *g) {
  data = d;
  cData = c;
  gData = g;

  u_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_data  = (double *)calloc(15 * data->numCells, sizeof(double));
  h_data    = (double *)calloc(data->numCells, sizeof(double));
  cMu_data  = (double *)calloc(46 * data->numCells, sizeof(double));
  gMu_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  cRho_data = (double *)calloc(46 * data->numCells, sizeof(double));
  gRho_data = (double *)calloc(21 * data->numCells, sizeof(double));

  u    = op_decl_dat(data->cells, 15, "double", u_data, "poisson_u");
  rhs  = op_decl_dat(data->cells, 15, "double", rhs_data, "poisson_rhs");
  h    = op_decl_dat(data->cells, 1, "double", h_data, "poisson_h");
  cMu  = op_decl_dat(data->cells, 46, "double", cMu_data, "poisson_cMu");
  gMu  = op_decl_dat(data->cells, 21, "double", gMu_data, "poisson_gMu");
  cRho = op_decl_dat(data->cells, 46, "double", cRho_data, "poisson_cRho");
  gRho = op_decl_dat(data->cells, 21, "double", gRho_data, "poisson_gRho");
}

ViscositySolve::~ViscositySolve() {
  free(u_data);
  free(rhs_data);
  free(h_data);
  free(cMu_data);
  free(gMu_data);
  free(cRho_data);
  free(gRho_data);

  destroy_vec(&b);
  destroy_vec(&x);
  KSPDestroy(&ksp);
  MatDestroy(&Amat);
}

void ViscositySolve::setDirichletBCs(int *d) {
  dirichlet = d;
}

void ViscositySolve::setNeumannBCs(int *n) {
  neumann = n;
}

void ViscositySolve::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double ViscositySolve::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}

void ViscositySolve::init() {
  create_vec(&b);
  create_vec(&x);
  create_shell_mat(&Amat);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  KSPSetOperators(ksp, Amat, Amat);
  KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 5e4);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

  op_par_loop(poisson_h, "poisson_h", data->cells,
              op_arg_dat(data->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(h, -1, OP_ID, 1, "double", OP_WRITE));
}

bool ViscositySolve::solve(op_dat b_dat, op_dat x_dat, double factor) {
  massFactor = factor;

  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_V), 15, data->nu, 0.0, cMu);
  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, data->nu, 0.0, gMu);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_V), 15, data->rho, 0.0, cRho);
  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, data->rho, 0.0, gRho);

  op_par_loop(viscosity_solve_apply_bc, "viscosity_solve_apply_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(gData->mD[0], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[1], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[2], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->sJ, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(h, 0, data->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(gData->tau, 0, data->bedge2cells, 3, "double", OP_READ),
              op_arg_dat(gMu, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(data->nu, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(gRho, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(bc_dat, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(b_dat, 0, data->bedge2cells, 15, "double", OP_INC));

  load_vec(&b, b_dat);
  load_vec(&x, x_dat);

  KSPSolve(ksp, b, x);
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

  return converged;
}

void ViscositySolve::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  copy_u(u_d);

  op_par_loop(viscosity_solve_0, "viscosity_solve_0", data->cells,
              op_arg_dat(cData->J,  -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->Dx, -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(cData->Dy, -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(cMu, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
              op_arg_gbl(&massFactor, 1, "double", OP_READ),
              op_arg_dat(cData->mm, -1, OP_ID, 15 * 15, "double", OP_READ),
              op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_WRITE));

  // gVP is stored in pDy (for the time being, will change later)
  op_par_loop(viscosity_solve_1, "viscosity_solve_1", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(gData->mD[0], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[0], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[1], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[1], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[2], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[2], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pD[0], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pD[0], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pD[1], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pD[1], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pD[2], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pD[2], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pDy[0], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pDy[0], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pDy[1], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pDy[1], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pDy[2], 0, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->pDy[2], 1, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->sJ, 0, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->sJ, 1, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(h, 0, data->edge2cells, 1, "double", OP_READ),
              op_arg_dat(h, 1, data->edge2cells, 1, "double", OP_READ),
              op_arg_dat(gData->tau, 0, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(gData->tau, 1, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(gMu, 0, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gMu, 1, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(data->nu, 0, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->nu, 1, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(gRho, 0, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gRho, 1, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(u, 0, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(u, 1, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(rhs, 0, data->edge2cells, 15, "double", OP_INC),
              op_arg_dat(rhs, 1, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(viscosity_solve_2, "viscosity_solve_2", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(gData->mD[0], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[1], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[2], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->sJ, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(h, 0, data->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(gData->tau, 0, data->bedge2cells, 3, "double", OP_READ),
              op_arg_dat(gMu, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(data->nu, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(gRho, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(u, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(rhs, 0, data->bedge2cells, 15, "double", OP_INC));

  copy_rhs(rhs_d);
}
