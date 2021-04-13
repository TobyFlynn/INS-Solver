#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "blas_calls.h"
#include "operators.h"

#include "kernels/gauss_tau.h"
#include "kernels/gauss_tau_bc.h"
#include "kernels/poisson_rhs_faces.h"
#include "kernels/poisson_rhs_bc.h"
#include "kernels/poisson_rhs_flux.h"
#include "kernels/poisson_rhs_J.h"

using namespace std;

Poisson_MF::Poisson_MF(INSData *nsData, CubatureData *cubData, GaussData *gaussData) : Poisson(nsData, cubData, gaussData) {
  createBCMatrix();

  u_data        = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_data      = (double *)calloc(15 * data->numCells, sizeof(double));
  tau_data      = (double *)calloc(3 * data->numCells, sizeof(double));
  gU_data       = (double *)calloc(21 * data->numCells, sizeof(double));
  uNumFlux_data = (double *)calloc(21 * data->numCells, sizeof(double));
  uFluxX_data   = (double *)calloc(21 * data->numCells, sizeof(double));
  uFluxY_data   = (double *)calloc(21 * data->numCells, sizeof(double));
  qx_data       = (double *)calloc(15 * data->numCells, sizeof(double));
  qy_data       = (double *)calloc(15 * data->numCells, sizeof(double));

  u        = op_decl_dat(data->cells, 15, "double", u_data, "poisson_u");
  rhs      = op_decl_dat(data->cells, 15, "double", rhs_data, "poisson_rhs");
  tau      = op_decl_dat(data->cells, 3, "double", tau_data, "poisson_tau");
  gU       = op_decl_dat(data->cells, 21, "double", gU_data, "poisson_gU");
  uNumFlux = op_decl_dat(data->cells, 21, "double", uNumFlux_data, "poisson_uNumFlux");
  uFluxX   = op_decl_dat(data->cells, 21, "double", uFluxX_data, "poisson_uFluxX");
  uFluxY   = op_decl_dat(data->cells, 21, "double", uFluxY_data, "poisson_uFluxY");
  qx       = op_decl_dat(data->cells, 15, "double", qx_data, "poisson_qx");
  qy       = op_decl_dat(data->cells, 15, "double", qy_data, "poisson_qy");

  // Calculate tau
  op_par_loop(gauss_tau, "gauss_tau", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->fscale, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(tau, -2, data->edge2cells, 3, "double", OP_INC));

  op_par_loop(gauss_tau_bc, "gauss_tau_bc", data->bedges,
              op_arg_dat(data->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->fscale, 0, data->bedge2cells, 15, "double", OP_READ),
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

  if(pBCMatInit)
    MatDestroy(&pBCMat);
}

bool Poisson_MF::solve(op_dat b_dat, op_dat x_dat, bool addMass, double factor) {
  massMat = addMass;
  massFactor = factor;
  Vec b;
  create_vec(&b);
  load_vec(&b, b_dat);

  Vec bc;
  create_vec(&bc, 21);
  load_vec(&bc, bc_dat, 21);

  // Calculate RHS for linear solve by applying the BCs
  Vec rhs;
  create_vec(&rhs);
  MatMultAdd(pBCMat, bc, b, rhs);

  Vec x;
  create_vec(&x);

  Mat Amat;
  create_shell_mat(&Amat);

  // Create PETSc Preconditioned Conjugate Gradient linear solver
  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetType(ksp, KSPCG);

  KSPSetOperators(ksp, Amat, Amat);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 1e4);
  // Solve
  KSPSolve(ksp, rhs, x);
  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  double residual;
  KSPGetResidualNorm(ksp, &residual);
  // Check that the solver converged
  bool converged = true;
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
  destroy_vec(&bc);
  destroy_vec(&rhs);
  MatDestroy(&Amat);

  return converged;
}

void Poisson_MF::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat (different depending on whether CPU or GPU)
  copy_u(u_d);

  gauss_interp_blas(data, u, gU);

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

  cub_grad(data, cData, u, qx, qy);

  poisson_rhs_blas1(data, this);

  op_par_loop(poisson_rhs_J, "poisson_rhs_J", data->cells,
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(qx, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(qy, -1, OP_ID, 15, "double", OP_RW));

  copy_rhs(rhs_d);
}

void Poisson_MF::createBCMatrix() {
  create_mat(&pBCMat, 15 * data->numCells, 21 * data->numCells, 15);
  pBCMatInit = true;
  double tol = 1e-15;

  double *gauss_sJ  = (double *)malloc(21 * op_get_size(data->cells) * sizeof(double));
  double *gauss_tau = (double *)malloc(3 * op_get_size(data->cells) * sizeof(double));
  double *gauss_mD[3];
  for(int i = 0; i < 3; i++) {
    gauss_mD[i]  = (double *)malloc(7 * 15 * op_get_size(data->cells) * sizeof(double));
    op_fetch_data(gData->mD[i], gauss_mD[i]);
  }
  op_fetch_data(gData->sJ, gauss_sJ);
  op_fetch_data(gData->tau, gauss_tau);

  // Create BCs matrix using Gauss data on boundary edges
  for(int i = 0; i < data->numBoundaryEdges; i++) {
    int element = data->bedge2cell_data[i];
    int bedgeType = data->bedge_type_data[i];
    int edge = data->bedgeNum_data[i];
    if(dirichlet[0] == bedgeType || dirichlet[1] == bedgeType || dirichlet[2] == bedgeType) {
      // Get data
      for(int j = 0; j < 7 * 15; j++) {
        int indT = (j % 7) * 15 + (j / 7);
        int col = element * 21 + edge * 7 + (j % 7);
        int row = element * 15 + (j / 7);
        double val;
        if(edge == 0) {
          val = gFInterp0[indT] * gaussW[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)] * gauss_tau[element * 3 + edge];
        } else if(edge == 1) {
          val = gFInterp1[indT] * gaussW[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)] * gauss_tau[element * 3 + edge];
        } else {
          val = gFInterp2[indT] * gaussW[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)] * gauss_tau[element * 3 + edge];
        }
        val -= gauss_mD[edge][element * 7 * 15 + indT] * gaussW[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)];
        if(abs(val) > tol)
          MatSetValues(pBCMat, 1, &row, 1, &col, &val, ADD_VALUES);
      }
    } else if(neumann[0] == bedgeType || neumann[1] == bedgeType || neumann[2] == bedgeType) {
      // Get data
      for(int j = 0; j < 7 * 15; j++) {
        int indT = (j % 7) * 15 + (j / 7);
        int col = element * 21 + edge * 7 + (j % 7);
        int row = element * 15 + (j / 7);
        double val;
        if(edge == 0) {
          val = gFInterp0[indT] * gaussW[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)];
        } else if(edge == 1) {
          val = gFInterp1[indT] * gaussW[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)];
        } else {
          val = gFInterp2[indT] * gaussW[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)];
        }
        if(abs(val) > tol)
          MatSetValues(pBCMat, 1, &row, 1, &col, &val, ADD_VALUES);
      }
    } else {
      cout << "UNDEFINED BOUNDARY EDGE" << endl;
      cout << "Element " << element << " Edge " << edge << " Type " << bedgeType << endl;
      cout << "D: " << dirichlet[0] << " " << dirichlet[1] << endl;
      cout << "N: " << neumann[0] << " " << neumann[1] << endl;
    }
  }

  free(gauss_sJ);
  free(gauss_tau);
  for(int i = 0; i < 3; i++) {
    free(gauss_mD[i]);
  }

  MatAssemblyBegin(pBCMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pBCMat, MAT_FINAL_ASSEMBLY);
}
