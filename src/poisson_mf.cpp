#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "blas_calls.h"
#include "operators.h"

#include "kernels/tau.h"
#include "kernels/tau_bc.h"
#include "kernels/poisson_rhs_faces.h"
#include "kernels/poisson_rhs_bc.h"
#include "kernels/poisson_rhs_flux.h"
#include "kernels/poisson_rhs_J.h"
#include "kernels/poisson_rhs_qbc.h"
#include "kernels/poisson_rhs_qflux.h"
#include "kernels/poisson_bc.h"
#include "kernels/poisson_bc_J.h"
#include "kernels/poisson_bc2.h"
#include "kernels/poisson_bc3.h"

using namespace std;

Poisson_MF::Poisson_MF(INSData *nsData, CubatureData *cubData, GaussData *gaussData) : Poisson(nsData, cubData, gaussData) {
  u_data         = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_data       = (double *)calloc(15 * data->numCells, sizeof(double));
  tau_data       = (double *)calloc(3 * data->numCells, sizeof(double));
  gU_data        = (double *)calloc(21 * data->numCells, sizeof(double));
  uNumFlux_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  uFluxX_data    = (double *)calloc(21 * data->numCells, sizeof(double));
  uFluxY_data    = (double *)calloc(21 * data->numCells, sizeof(double));
  dudx_data      = (double *)calloc(15 * data->numCells, sizeof(double));
  dudy_data      = (double *)calloc(15 * data->numCells, sizeof(double));
  qx_data        = (double *)calloc(15 * data->numCells, sizeof(double));
  qy_data        = (double *)calloc(15 * data->numCells, sizeof(double));
  gqx_data       = (double *)calloc(21 * data->numCells, sizeof(double));
  gqy_data       = (double *)calloc(21 * data->numCells, sizeof(double));
  qxNumFlux_data = (double *)calloc(21 * data->numCells, sizeof(double));
  qyNumFlux_data = (double *)calloc(21 * data->numCells, sizeof(double));
  qFlux_data     = (double *)calloc(21 * data->numCells, sizeof(double));
  gradx_data     = (double *)calloc(15 * data->numCells, sizeof(double));
  grady_data     = (double *)calloc(15 * data->numCells, sizeof(double));

  u         = op_decl_dat(data->cells, 15, "double", u_data, "poisson_u");
  rhs       = op_decl_dat(data->cells, 15, "double", rhs_data, "poisson_rhs");
  tau       = op_decl_dat(data->cells, 3, "double", tau_data, "poisson_tau");
  gU        = op_decl_dat(data->cells, 21, "double", gU_data, "poisson_gU");
  uNumFlux  = op_decl_dat(data->cells, 21, "double", uNumFlux_data, "poisson_uNumFlux");
  uFluxX    = op_decl_dat(data->cells, 21, "double", uFluxX_data, "poisson_uFluxX");
  uFluxY    = op_decl_dat(data->cells, 21, "double", uFluxY_data, "poisson_uFluxY");
  dudx      = op_decl_dat(data->cells, 15, "double", dudx_data, "poisson_dudx");
  dudy      = op_decl_dat(data->cells, 15, "double", dudy_data, "poisson_dudy");
  qx        = op_decl_dat(data->cells, 15, "double", qx_data, "poisson_qx");
  qy        = op_decl_dat(data->cells, 15, "double", qy_data, "poisson_qy");
  gqx       = op_decl_dat(data->cells, 21, "double", gqx_data, "poisson_gqx");
  gqy       = op_decl_dat(data->cells, 21, "double", gqy_data, "poisson_gqy");
  qxNumFlux = op_decl_dat(data->cells, 21, "double", qxNumFlux_data, "poisson_qxNumFlux");
  qyNumFlux = op_decl_dat(data->cells, 21, "double", qyNumFlux_data, "poisson_qyNumFlux");
  qFlux     = op_decl_dat(data->cells, 21, "double", qFlux_data, "poisson_qFlux");
  gradx     = op_decl_dat(data->cells, 15, "double", gradx_data, "poisson_gradx");
  grady     = op_decl_dat(data->cells, 15, "double", grady_data, "poisson_grady");

  // Calculate tau
  op_par_loop(tau, "tau", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->J, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->sJ, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(tau, -2, data->edge2cells, 3, "double", OP_INC));

  op_par_loop(tau_bc, "tau_bc", data->bedges,
              op_arg_dat(data->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->J, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->sJ, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(tau, 0, data->bedge2cells, 3, "double", OP_INC));
}

Poisson_MF::~Poisson_MF() {
  free(u_data);
  free(rhs_data);
  free(tau_data);
  free(gU_data);
  free(uNumFlux_data);
  free(uFluxX_data);
  free(uFluxY_data);
  free(dudx_data);
  free(dudy_data);
  free(qx_data);
  free(qy_data);
  free(gqx_data);
  free(gqy_data);
  free(qxNumFlux_data);
  free(qyNumFlux_data);
  free(qFlux_data);
}

bool Poisson_MF::solve(op_dat b_dat, op_dat x_dat, bool addMass, double factor) {
  massMat = addMass;
  massFactor = factor;

  applyBCs(b_dat);

  Vec b;
  create_vec(&b);
  load_vec(&b, b_dat);

  Vec x;
  create_vec(&x);

  Mat Amat;
  create_shell_mat(&Amat);

  // Create PETSc Preconditioned Conjugate Gradient linear solver
  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetType(ksp, KSPCG);
  // KSPSetType(ksp, KSPFGMRES);

  KSPSetOperators(ksp, Amat, Amat);
  KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 1e4);
  // Solve
  timer->startLinearSolveMF();
  KSPSolve(ksp, b, x);
  timer->endLinearSolveMF();
  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  double residual;
  KSPGetResidualNorm(ksp, &residual);
  // Check that the solver converged
  bool converged = true;
  cout << "Number of iterations for linear solver: " << numIt << endl;
  if(reason < 0) {
    converged = false;
    cout << "Number of iterations for linear solver: " << numIt << endl;
    cout << "Converged reason: " << reason << " Residual: " << residual << endl;
  }
  numberIter += numIt;
  solveCount++;

  // Get solution and free PETSc vectors and matrix
  Vec solution;
  KSPGetSolution(ksp, &solution);
  store_vec(&solution, x_dat);
  KSPDestroy(&ksp);
  destroy_vec(&b);
  destroy_vec(&x);
  MatDestroy(&Amat);

  return converged;
}

void Poisson_MF::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat (different depending on whether CPU or GPU)
  copy_u(u_d);
  timer->startLinearSolveMFRHS();
  gauss_interp_blas(data, u, gU);
  timer->endLinearSolveMFRHS();

  op_par_loop(poisson_rhs_faces, "poisson_rhs_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(gU, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(uNumFlux, -2, data->edge2cells, 21, "double", OP_INC));

  op_par_loop(poisson_rhs_bc, "poisson_rhs_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(gU, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(uNumFlux, 0, data->bedge2cells, 21, "double", OP_INC));

  op_par_loop(poisson_rhs_flux, "poisson_rhs_flux", data->cells,
              op_arg_dat(gData->nx, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(gData->ny, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(gData->sJ, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(uNumFlux, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(uFluxX, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(uFluxY, -1, OP_ID, 21, "double", OP_WRITE));

  timer->startLinearSolveMFRHS();
  cub_grad(data, cData, u, dudx, dudy);

  poisson_rhs_blas1(data, this);
  timer->endLinearSolveMFRHS();

  op_par_loop(poisson_rhs_J, "poisson_rhs_J", data->cells,
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(qx, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(qy, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(gradx, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(grady, -1, OP_ID, 15, "double", OP_RW));

  timer->startLinearSolveMFRHS();
  gauss_interp_blas(data, gradx, gqx);
  gauss_interp_blas(data, grady, gqy);
  timer->endLinearSolveMFRHS();

  op_par_loop(poisson_rhs_faces, "poisson_rhs_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(gqx, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(qxNumFlux, -2, data->edge2cells, 21, "double", OP_INC));

  op_par_loop(poisson_rhs_faces, "poisson_rhs_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(gqy, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(qyNumFlux, -2, data->edge2cells, 21, "double", OP_INC));

  op_par_loop(poisson_rhs_qbc, "poisson_rhs_qbc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&neumann[0], 1, "int", OP_READ),
              op_arg_gbl(&neumann[1], 1, "int", OP_READ),
              op_arg_gbl(&neumann[2], 1, "int", OP_READ),
              op_arg_dat(gqx, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(qxNumFlux, 0, data->bedge2cells, 21, "double", OP_INC));

  op_par_loop(poisson_rhs_qbc, "poisson_rhs_qbc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&neumann[0], 1, "int", OP_READ),
              op_arg_gbl(&neumann[1], 1, "int", OP_READ),
              op_arg_gbl(&neumann[2], 1, "int", OP_READ),
              op_arg_dat(gqy, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(qyNumFlux, 0, data->bedge2cells, 21, "double", OP_INC));

  op_par_loop(poisson_rhs_qflux, "poisson_rhs_qflux", data->cells,
              op_arg_dat(gData->nx, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(gData->ny, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(gData->sJ, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(tau, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(uNumFlux, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(qxNumFlux, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(qyNumFlux, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(qFlux, -1, OP_ID, 21, "double", OP_WRITE));

  timer->startLinearSolveMFRHS();
  cub_div(data, cData, qx, qy, rhs);

  poisson_rhs_blas2(data, this);
  timer->endLinearSolveMFRHS();

  if(massMat) {
    timer->startLinearSolveMFRHS();
    poisson_rhs_mass_blas(data, cData, this, massFactor);
    timer->endLinearSolveMFRHS();
  }

  copy_rhs(rhs_d);
}

void Poisson_MF::applyBCs(op_dat b_dat) {
  op_par_loop(poisson_bc, "poisson_bc", data->cells,
              op_arg_dat(gData->nx, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(gData->ny, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(gData->sJ, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(bc_dat, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(uFluxX, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(uFluxY, -1, OP_ID, 21, "double", OP_WRITE));

  poisson_bc_blas(data, this);

  op_par_loop(poisson_bc_J, "poisson_bc_J", data->cells,
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(gradx, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(grady, -1, OP_ID, 15, "double", OP_RW));

  cub_div(data, cData, gradx, grady, dudx);

  op_par_loop(poisson_bc2, "poisson_bc2", data->cells,
              op_arg_dat(gData->sJ, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(tau, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(bc_dat, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(uFluxX, -1, OP_ID, 21, "double", OP_WRITE));

  poisson_bc_blas2(data, this);

  op_par_loop(poisson_bc3, "poisson_bc3", data->cells,
              op_arg_dat(dudx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(b_dat, -1, OP_ID, 15, "double", OP_RW));
}
