#include "poisson.h"

#include <iostream>

#include "op_seq.h"
#include "blas_calls.h"
#include "operators.h"

#include "kernels/setup_poisson.h"
#include "kernels/set_tau.h"
#include "kernels/set_tau_bc.h"
#include "kernels/poisson_rhs_faces.h"
#include "kernels/poisson_rhs_bc.h"
#include "kernels/poisson_rhs_du.h"
#include "kernels/poisson_rhs_qbc.h"
#include "kernels/poisson_rhs_fluxq.h"
#include "kernels/poisson_rhs_J.h"

using namespace std;

PetscErrorCode matAMult(Mat A, Vec x, Vec y);

Poisson::Poisson(INSData *nsData) {
  data = nsData;
  // Allocate memory
  pTau_data      = (double *)malloc(15 * data->numCells * sizeof(double));
  pExRHS_data[0] = (double *)malloc(15 * data->numCells * sizeof(double));
  pExRHS_data[1] = (double *)malloc(15 * data->numCells * sizeof(double));
  pU_data        = (double *)malloc(15 * data->numCells * sizeof(double));
  pDu_data       = (double *)malloc(15 * data->numCells * sizeof(double));
  pFluxXu_data   = (double *)malloc(15 * data->numCells * sizeof(double));
  pFluxYu_data   = (double *)malloc(15 * data->numCells * sizeof(double));
  pDuDx_data     = (double *)malloc(15 * data->numCells * sizeof(double));
  pDuDy_data     = (double *)malloc(15 * data->numCells * sizeof(double));
  pFluxQ_data    = (double *)malloc(15 * data->numCells * sizeof(double));
  pDivQ_data     = (double *)malloc(15 * data->numCells * sizeof(double));
  pRHS_data      = (double *)malloc(15 * data->numCells * sizeof(double));
  // Declare OP2 dats
  pTau      = op_decl_dat(data->cells, 15, "double", pTau_data, "pTau");
  pExRHS[0] = op_decl_dat(data->cells, 15, "double", pExRHS_data[0], "pExRHS0");
  pExRHS[1] = op_decl_dat(data->cells, 15, "double", pExRHS_data[1], "pExRHS1");
  pU        = op_decl_dat(data->cells, 15, "double", pU_data, "pU");
  pDu       = op_decl_dat(data->cells, 15, "double", pDu_data, "pDu");
  pFluxXu   = op_decl_dat(data->cells, 15, "double", pFluxXu_data, "pFluxXu");
  pFluxYu   = op_decl_dat(data->cells, 15, "double", pFluxYu_data, "pFluxYu");
  pDuDx     = op_decl_dat(data->cells, 15, "double", pDuDx_data, "pDuDx");
  pDuDy     = op_decl_dat(data->cells, 15, "double", pDuDy_data, "pDuDy");
  pFluxQ    = op_decl_dat(data->cells, 15, "double", pFluxQ_data, "pFluxQ");
  pDivQ     = op_decl_dat(data->cells, 15, "double", pDivQ_data, "pDivQ");
  pRHS      = op_decl_dat(data->cells, 15, "double", pRHS_data, "pRHS");

  // Initialisation kernels
  op_par_loop(setup_poisson, "setup_poisson", data->cells,
              op_arg_dat(pTau, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(pExRHS[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(pExRHS[1], -1, OP_ID, 15, "double", OP_WRITE));

  op_par_loop(set_tau, "set_tau", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->J,  -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->sJ, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(pTau, -2, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(set_tau_bc, "set_tau_bc", data->bedges,
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->J,  0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->sJ, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(pTau, 0, data->bedge2cells, 15, "double", OP_INC));
}

Poisson::~Poisson() {
  // Free memory
  free(pTau_data);
  free(pExRHS_data[0]);
  free(pExRHS_data[1]);
  free(pU_data);
  free(pDu_data);
  free(pFluxXu_data);
  free(pFluxYu_data);
  free(pDuDx_data);
  free(pDuDy_data);
  free(pFluxQ_data);
  free(pDivQ_data);
}

void Poisson::solve(op_dat b_dat, op_dat x_dat) {
  Vec b;
  create_vec(&b);
  load_vec(&b, b_dat);

  Vec x;
  create_vec(&x);

  Mat Amat;
  MatCreateShell(PETSC_COMM_SELF, 15 * data->numCells, 15 * data->numCells, PETSC_DETERMINE, PETSC_DETERMINE, this, &Amat);
  MatShellSetOperation(Amat, MATOP_MULT, (void(*)(void))matAMult);

  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetType(ksp, KSPFGMRES);
  // KSPSetType(ksp, KSPCG);
  KSPSetOperators(ksp, Amat, Amat);

  // Solve
  KSPSolve(ksp, b, x);
  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  cout << "Number of iterations for linear solver: " << numIt << endl;
  cout << "Converged reason: " << reason << endl;

  Vec solution;
  KSPGetSolution(ksp, &solution);
  store_vec(&solution, x_dat);
  KSPDestroy(&ksp);
  MatDestroy(&Amat);
  destroy_vec(&b);
  destroy_vec(&x);
}

void Poisson::rhs(const double *u, double *rhs) {
  // Copy u to OP2 dat (different depending on whether CPU or GPU)
  copy_u(u);

  op_par_loop(poisson_rhs_faces, "poisson_rhs_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(pU, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(pExRHS[0], -2, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(poisson_rhs_bc, "poisson_rhs_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_dat(pU, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(pExRHS[0], 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(poisson_rhs_du, "poisson_rhs_du", data->cells,
              op_arg_dat(data->nx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ny, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->fscale, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(pU, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(pExRHS[0], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(pDu, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(pFluxXu, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(pFluxYu, -1, OP_ID, 15, "double", OP_WRITE));

  grad(data, pU, pDuDx, pDuDy);

  // qx and qy stored in pDuDx and pDuDy
  poisson_rhs_blas1(data, this);

  op_par_loop(poisson_rhs_faces, "poisson_rhs_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(pDuDx, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(pExRHS[0], -2, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(poisson_rhs_faces, "poisson_rhs_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(pDuDy, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(pExRHS[1], -2, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(poisson_rhs_qbc, "poisson_rhs_qbc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&neumann[0], 1, "int", OP_READ),
              op_arg_gbl(&neumann[1], 1, "int", OP_READ),
              op_arg_dat(pDuDx, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(pExRHS[0], 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(poisson_rhs_qbc, "poisson_rhs_qbc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&neumann[0], 1, "int", OP_READ),
              op_arg_gbl(&neumann[1], 1, "int", OP_READ),
              op_arg_dat(pDuDy, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(pExRHS[1], 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(poisson_rhs_fluxq, "poisson_rhs_fluxq", data->cells,
              op_arg_dat(data->nx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ny, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->fscale, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(pTau, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(pDu, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(pDuDx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(pDuDy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(pExRHS[0], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(pExRHS[1], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(pFluxQ, -1, OP_ID, 15, "double", OP_WRITE));

  div(data, pDuDx, pDuDy, pDivQ);

  poisson_rhs_blas2(data, this);

  op_par_loop(poisson_rhs_J, "poisson_rhs_J", data->cells,
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(pRHS, -1, OP_ID, 15, "double", OP_RW));

  // Different depending on whether CPU or GPU
  copy_rhs(rhs);
}

void Poisson::setDirichletBCs(int *d) {
  dirichlet = d;
}

void Poisson::setNeumannBCs(int *n) {
  neumann = n;
}