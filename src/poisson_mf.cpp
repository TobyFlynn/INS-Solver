#include "poisson.h"

#include "op_seq.h"

#include <iostream>

#include "blas_calls.h"
#include "operators.h"

using namespace std;

Poisson_MF::Poisson_MF(INSData *nsData, CubatureData *cubData, GaussData *gaussData) : Poisson(nsData, cubData, gaussData) {
  u_data      = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  gU_data     = (double *)calloc(21 * data->numCells, sizeof(double));
  gDudx_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  gDudy_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  fluxX_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  fluxY_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  flux_data   = (double *)calloc(21 * data->numCells, sizeof(double));
  dudx_data   = (double *)calloc(15 * data->numCells, sizeof(double));
  dudy_data   = (double *)calloc(15 * data->numCells, sizeof(double));
  qx_data     = (double *)calloc(15 * data->numCells, sizeof(double));
  qy_data     = (double *)calloc(15 * data->numCells, sizeof(double));
  tmp_u_data  = (double *)calloc(15 * data->numCells, sizeof(double));

  u      = op_decl_dat(data->cells, 15, "double", u_data, "poisson_u");
  rhs    = op_decl_dat(data->cells, 15, "double", rhs_data, "poisson_rhs");
  gU     = op_decl_dat(data->cells, 21, "double", gU_data, "poisson_gU");
  gDudx  = op_decl_dat(data->cells, 21, "double", gDudx_data, "poisson_gDudx");
  gDudy  = op_decl_dat(data->cells, 21, "double", gDudy_data, "poisson_gDudy");
  fluxX  = op_decl_dat(data->cells, 21, "double", fluxX_data, "poisson_fluxX");
  fluxY  = op_decl_dat(data->cells, 21, "double", fluxY_data, "poisson_fluxY");
  flux   = op_decl_dat(data->cells, 21, "double", flux_data, "poisson_flux");
  dudx   = op_decl_dat(data->cells, 15, "double", dudx_data, "poisson_dudx");
  dudy   = op_decl_dat(data->cells, 15, "double", dudy_data, "poisson_dudy");
  qx     = op_decl_dat(data->cells, 15, "double", qx_data, "poisson_qx");
  qy     = op_decl_dat(data->cells, 15, "double", qy_data, "poisson_qy");
  tmp_u  = op_decl_dat(data->cells, 15, "double", tmp_u_data, "poisson_tmp_u");
}

Poisson_MF::~Poisson_MF() {
  free(u_data);
  free(rhs_data);
  free(gU_data);
  free(gDudx_data);
  free(gDudy_data);
  free(fluxX_data);
  free(fluxY_data);
  free(flux_data);
  free(dudx_data);
  free(dudy_data);
  free(tmp_u_data);

  destroy_vec(&b);
  destroy_vec(&x);
  KSPDestroy(&ksp);
  MatDestroy(&Amat);
}

void Poisson_MF::init() {
  create_vec(&b);
  create_vec(&x);
  create_shell_mat(&Amat);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  KSPSetOperators(ksp, Amat, Amat);
  KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 1e5);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
}

bool Poisson_MF::solve(op_dat b_dat, op_dat x_dat, bool addMass, double factor) {
  massMat = addMass;
  massFactor = factor;
  scalarFactor = true;

  apply_bc(b_dat);

  load_vec(&b, b_dat);

  load_vec(&x, x_dat);

  // Solve
  timer->startLinearSolveMF();
  KSPSolve(ksp, b, x);
  timer->endLinearSolveMF();
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

  // Get solution and free PETSc vectors and matrix
  Vec solution;
  KSPGetSolution(ksp, &solution);
  store_vec(&solution, x_dat);

  return converged;
}

void Poisson_MF::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  copy_u(u_d);

  timer->startLinearSolveMFRHS();

  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, u, 0.0, gU);

  cub_grad(data, cData, u, dudx, dudy);
  inv_mass(data, dudx);
  inv_mass(data, dudy);

  if(massMat) {
    op_par_loop(poisson_mf_nu, "poisson_mf_nu", data->cells,
                op_arg_dat(data->nu, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(dudx, -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(dudy, -1, OP_ID, 15, "double", OP_RW));
  } else {
    op_par_loop(poisson_mf_rho, "poisson_mf_rho", data->cells,
                op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(dudx, -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(dudy, -1, OP_ID, 15, "double", OP_RW));
  }

  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, dudx, 0.0, gDudx);
  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, dudy, 0.0, gDudy);

  op_par_loop(poisson_mf_edges, "poisson_mf_edges", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(gData->sJ, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->nx, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->ny, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->tau, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(gU, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gDudx, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gDudy, -2, data->edge2cells, 21, "double", OP_READ),
              op_arg_dat(fluxX, -2, data->edge2cells, 21, "double", OP_INC),
              op_arg_dat(fluxY, -2, data->edge2cells, 21, "double", OP_INC),
              op_arg_dat(flux, -2, data->edge2cells, 21, "double", OP_INC));

  op_par_loop(poisson_mf_bedges, "poisson_mf_bedges", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(gData->sJ, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->nx, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->ny, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->tau, 0, data->bedge2cells, 3, "double", OP_READ),
              op_arg_dat(gU, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gDudx, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gDudy, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(fluxX, 0, data->bedge2cells, 21, "double", OP_INC),
              op_arg_dat(fluxY, 0, data->bedge2cells, 21, "double", OP_INC),
              op_arg_dat(flux, 0, data->bedge2cells, 21, "double", OP_INC));

  cub_grad_weak(data, cData, u, qx, qy);

  op2_gemv(false, 15, 21, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, fluxX, -1.0, qx);
  op2_gemv(false, 15, 21, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, fluxY, -1.0, qy);

  inv_mass(data, qx);
  inv_mass(data, qy);

  if(massMat) {
    op_par_loop(poisson_mf_nu, "poisson_mf_nu", data->cells,
                op_arg_dat(data->nu, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(dudx, -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(dudy, -1, OP_ID, 15, "double", OP_RW));
  } else {
    op_par_loop(poisson_mf_rho, "poisson_mf_rho", data->cells,
                op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(dudx, -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(dudy, -1, OP_ID, 15, "double", OP_RW));
  }

  cub_div_weak(data, cData, qx, qy, rhs);

  op2_gemv(false, 15, 21, -1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, flux, 1.0, rhs);

  op_par_loop(poisson_mf_zero, "poisson_mf_zero", data->cells,
              op_arg_dat(fluxX, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(fluxY, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(flux, -1, OP_ID, 21, "double", OP_WRITE));

  if(massMat) {
    op_par_loop(poisson_mf_mm_rho, "poisson_mf_mm_rho", data->cells,
                op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(tmp_u, -1, OP_ID, 15, "double", OP_WRITE));

    op2_gemv_batch(false, 15, 15, massFactor, cData->mm, 15, tmp_u, 1.0, rhs);
  }

  timer->endLinearSolveMFRHS();

  copy_rhs(rhs_d);
}

void Poisson_MF::apply_bc(op_dat b) {
  op_par_loop(poisson_mf_bc0, "poisson_mf_bc0", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(gData->sJ, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->nx, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->ny, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->tau, 0, data->bedge2cells, 3, "double", OP_READ),
              op_arg_dat(bc_dat, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(fluxX, 0, data->bedge2cells, 21, "double", OP_INC),
              op_arg_dat(fluxY, 0, data->bedge2cells, 21, "double", OP_INC),
              op_arg_dat(flux, 0, data->bedge2cells, 21, "double", OP_INC));

  op2_gemv(false, 15, 21, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, fluxX, 0.0, qx);
  op2_gemv(false, 15, 21, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, fluxY, 0.0, qy);

  if(massMat) {
    op_par_loop(poisson_mf_nu, "poisson_mf_nu", data->cells,
                op_arg_dat(data->nu, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(qx, -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(qy, -1, OP_ID, 15, "double", OP_RW));
  } else {
    op_par_loop(poisson_mf_rho, "poisson_mf_rho", data->cells,
                op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(qx, -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(qy, -1, OP_ID, 15, "double", OP_RW));
  }

  cub_div_weak(data, cData, qx, qy, rhs);

  op2_gemv(false, 15, 21, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, flux, 1.0, rhs);

  op_par_loop(poisson_mf_bc1, "poisson_mf_bc1", data->cells,
              op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(b, -1, OP_ID, 15, "double", OP_RW));

  op_par_loop(poisson_mf_zero, "poisson_mf_zero", data->cells,
              op_arg_dat(fluxX, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(fluxY, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(flux, -1, OP_ID, 21, "double", OP_WRITE));
}
