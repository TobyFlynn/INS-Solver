#include "poisson.h"

#include <iostream>
#include <unistd.h>

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

Poisson::Poisson(INSData *nsData, CubatureData *cubData, GaussData *gaussData) {
  data = nsData;
  cData = cubData;
  gData = gaussData;
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

  Vec bc;
  create_vec(&bc, 21);
  load_vec(&bc, dBC, 21);

  Vec rhs;
  create_vec(&rhs);
  MatMultAdd(pBCMat, bc, b, rhs);

  Vec x;
  create_vec(&x);

  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  // if(method) {
    // KSPSetType(ksp, KSPFGMRES);
  // } else {
    KSPSetType(ksp, KSPCG);
  // }
  // KSPSetType(ksp, KSPGMRES);
  // KSPSetType(ksp, KSPCG);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCICC);

  Mat op;
  MatCreate(PETSC_COMM_SELF, &op);
  MatSetSizes(op, PETSC_DECIDE, PETSC_DECIDE, 15 * data->numCells, 15 * data->numCells);
  MatSetType(op, MATSEQAIJ);
  MatSeqAIJSetPreallocation(op, 15 * 4, NULL);
  MatAssemblyBegin(op, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(op, MAT_FINAL_ASSEMBLY);
  if(addMass) {
    MatCopy(pMat, op, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(op, factor, pMMat, DIFFERENT_NONZERO_PATTERN);
  } else {
    MatCopy(pMat, op, DIFFERENT_NONZERO_PATTERN);
  }

  KSPSetOperators(ksp, op, op);
  KSPSetTolerances(ksp, 1e-12, 1e-50, 1e5, 1e5);
  // Solve
  KSPSolve(ksp, rhs, x);
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
  MatDestroy(&op);
}

void Poisson::setDirichletBCs(int *d) {
  dirichlet = d;
}

void Poisson::setNeumannBCs(int *n) {
  neumann = n;
}

void Poisson::setBCValues(op_dat d_dat) {
  dBC = d_dat;
}

void Poisson::createMatrix() {
  MatCreate(PETSC_COMM_SELF, &pMat);
  MatSetSizes(pMat, PETSC_DECIDE, PETSC_DECIDE, 15 * data->numCells, 15 * data->numCells);
  MatSetType(pMat, MATSEQAIJ);
  MatSeqAIJSetPreallocation(pMat, 15 * 4, NULL);
  double tol = 1e-15;
  // double tol = 0.0;
  // Add cubature OP
  double *cub_OP = (double *)malloc(15 * 15 * op_get_size(data->cells) * sizeof(double));
  op_fetch_data(cData->OP, cub_OP);
  for(int i = 0; i < data->numCells; i++) {
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = i * 15 + m;
        int col = i * 15 + n;
        int colInd = n * 15 + m;
        double val = cub_OP[i * 15 * 15 + colInd];
        if(abs(val) > tol)
          MatSetValues(pMat, 1, &row, 1, &col, &val, ADD_VALUES);
      }
    }
  }
  free(cub_OP);

  double *gauss_OP[3];
  double *gauss_OPf[3];

  for(int i = 0; i < 3; i++) {
    gauss_OP[i] = (double *)malloc(15 * 15 * op_get_size(data->cells) * sizeof(double));
    gauss_OPf[i] = (double *)malloc(15 * 15 * op_get_size(data->cells) * sizeof(double));
    op_fetch_data(gData->OP[i], gauss_OP[i]);
    op_fetch_data(gData->OPf[i], gauss_OPf[i]);
  }

  // Gauss OP and OPf
  for(int i = 0; i < data->numEdges; i++) {
    int leftElement = data->edge2cell_data[i * 2];
    int rightElement = data->edge2cell_data[i * 2 + 1];
    int leftEdge = data->edgeNum_data[i * 2];
    int rightEdge = data->edgeNum_data[i * 2 + 1];
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = leftElement * 15 + m;
        int col = leftElement * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gauss_OP[leftEdge][leftElement * 15 * 15 + colInd];
        if(abs(val) > tol)
          MatSetValues(pMat, 1, &row, 1, &col, &val, ADD_VALUES);
      }
    }
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = rightElement * 15 + m;
        int col = rightElement * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gauss_OP[rightEdge][rightElement * 15 * 15 + colInd];
        if(abs(val) > tol)
          MatSetValues(pMat, 1, &row, 1, &col, &val, ADD_VALUES);
      }
    }

    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = leftElement * 15 + m;
        int col = rightElement * 15 + n;
        int colInd = n * 15 + m;
        double val = -0.5 * gauss_OPf[leftEdge][leftElement * 15 * 15 + colInd];
        if(abs(val) > tol)
          MatSetValues(pMat, 1, &row, 1, &col, &val, ADD_VALUES);
      }
    }

    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = rightElement * 15 + m;
        int col = leftElement * 15 + n;
        int colInd = n * 15 + m;
        double val = -0.5 * gauss_OPf[rightEdge][rightElement * 15 * 15 + colInd];
        if(abs(val) > tol)
          MatSetValues(pMat, 1, &row, 1, &col, &val, ADD_VALUES);
      }
    }
  }

  // Gauss on boundary OP
  for(int i = 0; i < data->numBoundaryEdges; i++) {
    int element = data->bedge2cell_data[i];
    int bedgeType = data->bedge_type_data[i];
    int edge = data->bedgeNum_data[i];
    if(dirichlet[0] == bedgeType || dirichlet[1] == bedgeType) {
      // Convert data to row major format
      for(int m = 0; m < 15; m++) {
        for(int n = 0; n < 15; n++) {
          int row = element * 15 + m;
          int col = element * 15 + n;
          int colInd = n * 15 + m;
          double val = gauss_OP[edge][element * 15 * 15 + colInd];
          if(abs(val) > tol)
            MatSetValues(pMat, 1, &row, 1, &col, &val, ADD_VALUES);
        }
      }
    }
  }

  for(int i = 0; i < 3; i++) {
    free(gauss_OP[i]);
    free(gauss_OPf[i]);
  }

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
  // PetscViewer pv;
  // PetscViewerDrawOpen(PETSC_COMM_SELF, NULL, NULL, PETSC_DECIDE, PETSC_DECIDE, 500, 500, &pv);
  // MatView(pMat, pv);
}

void Poisson::createMassMatrix() {
  MatCreate(PETSC_COMM_SELF, &pMMat);
  MatSetSizes(pMMat, PETSC_DECIDE, PETSC_DECIDE, 15 * data->numCells, 15 * data->numCells);
  MatSetType(pMMat, MATSEQAIJ);
  MatSeqAIJSetPreallocation(pMMat, 15, NULL);

  // Add cubature OP
  double *cub_MM = (double *)malloc(15 * 15 * op_get_size(data->cells) * sizeof(double));
  op_fetch_data(cData->mm, cub_MM);
  for(int i = 0; i < data->numCells; i++) {
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = i * 15 + m;
        int col = i * 15 + n;
        int colInd = n * 15 + m;
        double val = cub_MM[i * 15 * 15 + colInd];
        MatSetValues(pMMat, 1, &row, 1, &col, &val, ADD_VALUES);
      }
    }
  }
  free(cub_MM);

  MatAssemblyBegin(pMMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMMat, MAT_FINAL_ASSEMBLY);
}

void Poisson::createBCMatrix() {
  MatCreate(PETSC_COMM_SELF, &pBCMat);
  MatSetSizes(pBCMat, PETSC_DECIDE, PETSC_DECIDE, 15 * data->numCells, 21 * data->numCells);
  MatSetUp(pBCMat);

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

  for(int i = 0; i < data->numBoundaryEdges; i++) {
    int element = data->bedge2cell_data[i];
    int bedgeType = data->bedge_type_data[i];
    int edge = data->bedgeNum_data[i];
    if(dirichlet[0] == bedgeType || dirichlet[1] == bedgeType) {
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
    } else if(neumann[0] == bedgeType || neumann[1] == bedgeType) {
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
