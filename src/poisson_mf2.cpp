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

Poisson_MF2::Poisson_MF2(INSData *nsData, CubatureData *cubData, GaussData *gaussData) : Poisson(nsData, cubData, gaussData) {
  u_data      = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  op1_data    = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  op2_data[0] = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  op2_data[1] = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  op2_data[2] = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  u_t_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_t_data  = (double *)calloc(15 * data->numCells, sizeof(double));

  u     = op_decl_dat(data->cells, 15, "double", u_data, "poisson_u");
  rhs   = op_decl_dat(data->cells, 15, "double", rhs_data, "poisson_rhs");
  u_t   = op_decl_dat(data->cells, 15, "double", u_t_data, "poisson_u_t");
  rhs_t = op_decl_dat(data->cells, 15, "double", rhs_t_data, "poisson_rhs_t");
}

Poisson_MF2::~Poisson_MF2() {
  free(u_data);
  free(rhs_data);
  free(op1_data);
  free(op2_data[0]);
  free(op2_data[1]);
  free(op2_data[2]);
  free(u_t_data);
  free(rhs_t_data);
}

bool Poisson_MF2::solve(op_dat b_dat, op_dat x_dat, bool addMass, double factor) {
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
  // KSPSetType(ksp, KSPFGMRES);

  KSPSetOperators(ksp, Amat, Amat);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 1e5);
  // Solve
  timer->startLinearSolveMF();
  KSPSolve(ksp, rhs, x);
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
  destroy_vec(&rhs);
  destroy_vec(&bc);
  MatDestroy(&Amat);

  return converged;
}

void Poisson_MF2::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat (different depending on whether CPU or GPU)
  copy_u(u_d);

  poisson_mf2_blas(data, this);

  copy_rhs(rhs_d);
}

void Poisson_MF2::setOp() {
  double tol = 1e-15;

  // Add cubature OP to Poisson matrix
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
        if(abs(val) > tol) {
          int ind = i * 15 * 15 + m * 15 + n;
          op1_data[ind] += val;
        }
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

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < data->numEdges; i++) {
    int leftElement = data->edge2cell_data[i * 2];
    int rightElement = data->edge2cell_data[i * 2 + 1];
    int leftEdge = data->edgeNum_data[i * 2];
    int rightEdge = data->edgeNum_data[i * 2 + 1];
    // Gauss OP
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = leftElement * 15 + m;
        int col = leftElement * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gauss_OP[leftEdge][leftElement * 15 * 15 + colInd];
        if(abs(val) > tol) {
          int ind = leftElement * 15 * 15 + m * 15 + n;
          op1_data[ind] += val;
        }
      }
    }
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = rightElement * 15 + m;
        int col = rightElement * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gauss_OP[rightEdge][rightElement * 15 * 15 + colInd];
        if(abs(val) > tol) {
          int ind = rightElement * 15 * 15 + m * 15 + n;
          op1_data[ind] += val;
        }
      }
    }

    // Gauss OPf
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = leftElement * 15 + m;
        int col = rightElement * 15 + n;
        int colInd = n * 15 + m;
        double val = -0.5 * gauss_OPf[leftEdge][leftElement * 15 * 15 + colInd];
        if(abs(val) > tol) {
          int ind = leftElement * 15 * 15 + m * 15 + n;
          op2_data[leftEdge][ind] += val;
        }
      }
    }
    // Convert data to row major format
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int row = rightElement * 15 + m;
        int col = leftElement * 15 + n;
        int colInd = n * 15 + m;
        double val = -0.5 * gauss_OPf[rightEdge][rightElement * 15 * 15 + colInd];
        if(abs(val) > tol) {
          int ind = rightElement * 15 * 15 + m * 15 + n;
          op2_data[rightEdge][ind] += val;
        }
      }
    }
  }

  // Add Gauss OP for boundary edges
  for(int i = 0; i < data->numBoundaryEdges; i++) {
    int element = data->bedge2cell_data[i];
    int bedgeType = data->bedge_type_data[i];
    int edge = data->bedgeNum_data[i];
    if(dirichlet[0] == bedgeType || dirichlet[1] == bedgeType || dirichlet[2] == bedgeType) {
      // Convert data to row major format
      for(int m = 0; m < 15; m++) {
        for(int n = 0; n < 15; n++) {
          int row = element * 15 + m;
          int col = element * 15 + n;
          int colInd = n * 15 + m;
          double val = gauss_OP[edge][element * 15 * 15 + colInd];
          if(abs(val) > tol) {
            int ind = element * 15 * 15 + m * 15 + n;
            op1_data[ind] += val;
          }
        }
      }
    }
  }

  for(int i = 0; i < 3; i++) {
    free(gauss_OP[i]);
    free(gauss_OPf[i]);
  }

  op1    = op_decl_dat(data->cells, 15 * 15, "double", op1_data, "poisson_op1");
  op2[0] = op_decl_dat(data->cells, 15 * 15, "double", op2_data[0], "poisson_op20");
  op2[1] = op_decl_dat(data->cells, 15 * 15, "double", op2_data[1], "poisson_op21");
  op2[2] = op_decl_dat(data->cells, 15 * 15, "double", op2_data[2], "poisson_op22");

  // MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  // MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
  // PetscViewer pv;
  // PetscViewerDrawOpen(PETSC_COMM_SELF, NULL, NULL, PETSC_DECIDE, PETSC_DECIDE, 500, 500, &pv);
  // MatView(pMat, pv);
}

void Poisson_MF2::createBCMatrix() {
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
          val = gFInterp0_g[indT] * gaussW_g[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)] * gauss_tau[element * 3 + edge];
        } else if(edge == 1) {
          val = gFInterp1_g[indT] * gaussW_g[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)] * gauss_tau[element * 3 + edge];
        } else {
          val = gFInterp2_g[indT] * gaussW_g[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)] * gauss_tau[element * 3 + edge];
        }
        val -= gauss_mD[edge][element * 7 * 15 + indT] * gaussW_g[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)];
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
          val = gFInterp0_g[indT] * gaussW_g[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)];
        } else if(edge == 1) {
          val = gFInterp1_g[indT] * gaussW_g[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)];
        } else {
          val = gFInterp2_g[indT] * gaussW_g[j % 7] * gauss_sJ[element * 21 + edge * 7 + (j % 7)];
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
