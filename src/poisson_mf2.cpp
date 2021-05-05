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
#include "kernels/poisson_mf2.h"
#include "kernels/poisson_mf2_mass.h"
#include "kernels/poisson_mf2_faces.h"
#include "kernels/poisson_mf2_op.h"
#include "kernels/poisson_mf2_opf.h"
#include "kernels/poisson_mf2_opbf.h"
#include "kernels/poisson_mf2_bc.h"
#include "kernels/poisson_mf2_apply_bc.h"

using namespace std;

Poisson_MF2::Poisson_MF2(INSData *nsData, CubatureData *cubData, GaussData *gaussData) : Poisson(nsData, cubData, gaussData) {
  use_blas = false;

  u_data      = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  op1_data    = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  op2_data[0] = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  op2_data[1] = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  op2_data[2] = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  op_bc_data  = (double *)calloc(7 * 15 * data->numBoundaryEdges, sizeof(double));
  u_t_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_t_data  = (double *)calloc(15 * data->numCells, sizeof(double));

  u      = op_decl_dat(data->cells, 15, "double", u_data, "poisson_u");
  rhs    = op_decl_dat(data->cells, 15, "double", rhs_data, "poisson_rhs");
  u_t    = op_decl_dat(data->cells, 15, "double", u_t_data, "poisson_u_t");
  rhs_t  = op_decl_dat(data->cells, 15, "double", rhs_t_data, "poisson_rhs_t");
  op1    = op_decl_dat(data->cells, 15 * 15, "double", op1_data, "poisson_op1");
  op2[0] = op_decl_dat(data->cells, 15 * 15, "double", op2_data[0], "poisson_op20");
  op2[1] = op_decl_dat(data->cells, 15 * 15, "double", op2_data[1], "poisson_op21");
  op2[2] = op_decl_dat(data->cells, 15 * 15, "double", op2_data[2], "poisson_op22");
  op_bc  = op_decl_dat(data->bedges, 7 * 15, "double", op_bc_data, "poisson_op_bc");
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
  free(op_bc_data);
}

bool Poisson_MF2::solve(op_dat b_dat, op_dat x_dat, bool addMass, double factor) {
  massMat = addMass;
  massFactor = factor;

  op_par_loop(poisson_mf2_apply_bc, "poisson_mf2_apply_bc", data->bedges,
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(bc_dat, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(b_dat, 0, data->bedge2cells, 15, "double", OP_INC));

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
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 1e4);
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

void Poisson_MF2::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat (different depending on whether CPU or GPU)
  copy_u(u_d);

  if(use_blas) {
    // op2_gemv_batch(true, 15, 15, 1.0, op1, 15, u, 0.0, rhs);
    // if(massMat) {
    //   op2_gemv_batch(false, 15, 15, massFactor, cData->mm, 15, u, 1.0, rhs);
    // }
    // poisson_mf2_blas(data, this, cData, massMat, massFactor);
  } else {
    if(massMat) {
      op_par_loop(poisson_mf2_mass, "poisson_mf2_mass", data->cells,
                  op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_READ),
                  op_arg_gbl(&massFactor, 1, "double", OP_READ),
                  op_arg_dat(cData->mm, -1, OP_ID, 15 * 15, "double", OP_READ),
                  op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_WRITE));
    } else {
      op_par_loop(poisson_mf2, "poisson_mf2", data->cells,
                  op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_READ),
                  op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_WRITE));
    }
    op_par_loop(poisson_mf2_faces, "poisson_mf2_faces", data->edges,
                op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
                op_arg_dat(u, 0, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(op2[0], 0, data->edge2cells, 15 * 15, "double", OP_READ),
                op_arg_dat(op2[1], 0, data->edge2cells, 15 * 15, "double", OP_READ),
                op_arg_dat(op2[2], 0, data->edge2cells, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, 0, data->edge2cells, 15, "double", OP_INC),
                op_arg_dat(u, 1, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(op2[0], 1, data->edge2cells, 15 * 15, "double", OP_READ),
                op_arg_dat(op2[1], 1, data->edge2cells, 15 * 15, "double", OP_READ),
                op_arg_dat(op2[2], 1, data->edge2cells, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, 1, data->edge2cells, 15, "double", OP_INC));
  }

  copy_rhs(rhs_d);
}

void Poisson_MF2::setOp() {
  double tol = 1e-15;

  op_par_loop(poisson_mf2_op, "poisson_mf2_op", data->cells,
              op_arg_dat(cData->OP, -1, OP_ID, 15 * 15, "double", OP_READ),
              op_arg_gbl(&tol, 1, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_WRITE));

  op_par_loop(poisson_mf2_opf, "poisson_mf2_opf", data->edges,
              op_arg_gbl(&tol, 1, "double", OP_READ),
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(gData->OP[0], 0, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OP[1], 0, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OP[2], 0, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[0], 0, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[1], 0, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[2], 0, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(op2[0], 0, data->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(op2[1], 0, data->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(op2[2], 0, data->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(op1, 0, data->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(gData->OP[0], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OP[1], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OP[2], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[0], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[1], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[2], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(op2[0], 1, data->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(op2[1], 1, data->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(op2[2], 1, data->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(op1, 1, data->edge2cells, 15 * 15, "double", OP_INC));

  op_par_loop(poisson_mf2_opbf, "poisson_mf2_opbf", data->bedges,
              op_arg_gbl(&tol, 1, "double", OP_READ),
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(gData->OP[0], 0, data->bedge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OP[1], 0, data->bedge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OP[2], 0, data->bedge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(op1, 0, data->bedge2cells, 15 * 15, "double", OP_INC));
}

void Poisson_MF2::setBCOP() {
  double tol = 1e-15;
  // If not dirichlet BC, kernel will assume it is a neumann bc
  op_par_loop(poisson_mf2_bc, "poisson_mf2_bc", data->bedges,
              op_arg_gbl(&tol, 1, "double", OP_READ),
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(gData->mD[0], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[1], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->mD[2], 0, data->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(gData->sJ, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gData->tau, 0, data->bedge2cells, 3, "double", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, 7 * 15, "double", OP_INC));
}
