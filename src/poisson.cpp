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

void Poisson::solve(op_dat b_dat, op_dat x_dat, bool method, bool addMass, double factor) {
  // op_fetch_data_hdf5_file(dBC, "p.h5");
  massMat = addMass;
  massFactor = factor;
  Vec b;
  create_vec(&b);
  load_vec(&b, b_dat);

  Vec x;
  create_vec(&x);

  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  // if(method) {
    KSPSetType(ksp, KSPFGMRES);
  // } else {
  //   KSPSetType(ksp, KSPCG);
  // }
  // KSPSetType(ksp, KSPCG);
  // PC pc;
  // KSPGetPC(ksp, &pc);
  // PCSetType(pc, PCICC);
  // KSPSetPC(ksp, pc);
  // KSPSetPCSide(ksp, PC_RIGHT);
  KSPSetOperators(ksp, pMat, pMat);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 1e4);

  // Solve
  KSPSolve(ksp, b, x);
  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  double residual;
  KSPGetResidualNorm(ksp, &residual);
  cout << "Number of iterations for linear solver: " << numIt << endl;
  cout << "Converged reason: " << reason << " Residual: " << residual << endl;

  Vec solution;
  KSPGetSolution(ksp, &solution);
  store_vec(&solution, x_dat);
  KSPDestroy(&ksp);
  destroy_vec(&b);
  destroy_vec(&x);
}

void Poisson::setDirichletBCs(int *d, op_dat d_dat) {
  dirichlet = d;
  dBC = d_dat;
}

void Poisson::createMatrix() {
  MatCreate(PETSC_COMM_SELF, &pMat);
  MatSetSizes(pMat, PETSC_DECIDE, PETSC_DECIDE, 15 * data->numCells, 15 * data->numCells);
  MatSetUp(pMat);

  // TODO: Insert elements

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
}

void Poisson::createBCMatrix() {
  MatCreate(PETSC_COMM_SELF, &pBCMat);
  MatSetSizes(pBCMat, PETSC_DECIDE, PETSC_DECIDE, 15 * data->numCells, 15 * data->numCells);
  MatSetUp(pBCMat);

  // TODO: Insert elements

  MatAssemblyBegin(pBCMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pBCMat, MAT_FINAL_ASSEMBLY);
}
