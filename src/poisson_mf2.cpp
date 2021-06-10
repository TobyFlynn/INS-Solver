#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "blas_calls.h"
#include "operators.h"

using namespace std;

Poisson_MF2::Poisson_MF2(INSData *nsData, CubatureData *cubData, GaussData *gaussData) : Poisson(nsData, cubData, gaussData) {
  u_data      = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  op1_data    = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  op2_data[0] = (double *)calloc(15 * 15 * data->numEdges, sizeof(double));
  op2_data[1] = (double *)calloc(15 * 15 * data->numEdges, sizeof(double));
  op_bc_data  = (double *)calloc(7 * 15 * data->numBoundaryEdges, sizeof(double));
  u_t_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_t_data  = (double *)calloc(15 * data->numCells, sizeof(double));

  u      = op_decl_dat(data->cells, 15, "double", u_data, "poisson_u");
  rhs    = op_decl_dat(data->cells, 15, "double", rhs_data, "poisson_rhs");
  u_t    = op_decl_dat(data->cells, 15, "double", u_t_data, "poisson_u_t");
  rhs_t  = op_decl_dat(data->cells, 15, "double", rhs_t_data, "poisson_rhs_t");
  op1    = op_decl_dat(data->cells, 15 * 15, "double", op1_data, "poisson_op1");
  op2[0] = op_decl_dat(data->edges, 15 * 15, "double", op2_data[0], "poisson_op20");
  op2[1] = op_decl_dat(data->edges, 15 * 15, "double", op2_data[1], "poisson_op21");
  op_bc  = op_decl_dat(data->bedges, 7 * 15, "double", op_bc_data, "poisson_op_bc");
}

Poisson_MF2::~Poisson_MF2() {
  free(u_data);
  free(rhs_data);
  free(op1_data);
  free(op2_data[0]);
  free(op2_data[1]);
  free(u_t_data);
  free(rhs_t_data);
  free(op_bc_data);

  destroy_vec(&b);
  destroy_vec(&x);
  KSPDestroy(&ksp);
  MatDestroy(&Amat);
}

void Poisson_MF2::init() {
  setOp();
  setBCOP();

  create_vec(&b);
  create_vec(&x);
  create_shell_mat(&Amat);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  KSPSetOperators(ksp, Amat, Amat);
  KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 2e4);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
}

bool Poisson_MF2::solve(op_dat b_dat, op_dat x_dat, bool addMass, double factor) {
  massMat = addMass;
  massFactor = factor;
  scalarFactor = true;

  if(massMat) {
    // Viscosity Linear solve
    op_par_loop(poisson_mf2_apply_bc_vis, "poisson_mf2_apply_bc_vis", data->bedges,
                op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(op_bc, -1, OP_ID, 7 * 15, "double", OP_READ),
                op_arg_dat(data->nu, 0, data->bedge2cells, 15, "double", OP_READ),
                op_arg_dat(data->rho, 0, data->bedge2cells, 15, "double", OP_READ),
                op_arg_dat(bc_dat, 0, data->bedge2cells, 21, "double", OP_READ),
                op_arg_dat(b_dat, 0, data->bedge2cells, 15, "double", OP_INC));
  } else {
    // Pressure Linear solve
    op_par_loop(poisson_mf2_apply_bc, "poisson_mf2_apply_bc", data->bedges,
                op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(op_bc, -1, OP_ID, 7 * 15, "double", OP_READ),
                op_arg_dat(bc_dat, 0, data->bedge2cells, 21, "double", OP_READ),
                op_arg_dat(b_dat, 0, data->bedge2cells, 15, "double", OP_INC));
  }

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

void Poisson_MF2::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  copy_u(u_d);

  timer->startLinearSolveMFRHS();

  if(massMat) {
    op_par_loop(poisson_mf2_mass, "poisson_mf2_mass", data->cells,
                op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
                op_arg_gbl(&massFactor, 1, "double", OP_READ),
                op_arg_dat(cData->mm, -1, OP_ID, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_WRITE));

    op_par_loop(poisson_mf2_vis, "poisson_mf2_vis", data->cells,
                op_arg_dat(data->nu, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_RW));

    op_par_loop(poisson_mf2_faces_vis, "poisson_mf2_faces_vis", data->edges,
                op_arg_dat(data->nu, 0, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(data->rho, 0, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(u, 0, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(op2[0], -1, OP_ID, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, 0, data->edge2cells, 15, "double", OP_INC),
                op_arg_dat(data->nu, 1, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(data->rho, 1, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(u, 1, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(op2[1], -1, OP_ID, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, 1, data->edge2cells, 15, "double", OP_INC));
  } else {
    op_par_loop(poisson_mf2, "poisson_mf2", data->cells,
                op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_WRITE));

    op_par_loop(poisson_mf2_faces, "poisson_mf2_faces", data->edges,
                op_arg_dat(u, 0, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(op2[0], -1, OP_ID, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, 0, data->edge2cells, 15, "double", OP_INC),
                op_arg_dat(u, 1, data->edge2cells, 15, "double", OP_READ),
                op_arg_dat(op2[1], -1, OP_ID, 15 * 15, "double", OP_READ),
                op_arg_dat(rhs, 1, data->edge2cells, 15, "double", OP_INC));
  }

  timer->endLinearSolveMFRHS();

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
              op_arg_dat(op2[0], -1, OP_ID, 15 * 15, "double", OP_INC),
              op_arg_dat(op1, 0, data->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(gData->OP[0], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OP[1], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OP[2], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[0], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[1], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(gData->OPf[2], 1, data->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, 15 * 15, "double", OP_INC),
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
