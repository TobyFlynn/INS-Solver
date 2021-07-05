#include "poisson.h"

#include "op_seq.h"

#include <iostream>

#include "dg_constants.h"
#include "dg_blas_calls.h"

extern DGConstants *constants;

using namespace std;

PoissonSolve::PoissonSolve(DGMesh *m, INSData *d, bool p) {
  mesh = m;
  data = d;
  precondition = p;

  u_data      = (double *)calloc(15 * mesh->numCells, sizeof(double));
  rhs_data    = (double *)calloc(15 * mesh->numCells, sizeof(double));
  h_data      = (double *)calloc(mesh->numCells, sizeof(double));
  op1_data    = (double *)calloc(15 * 15 * mesh->numCells, sizeof(double));
  op2_data[0] = (double *)calloc(15 * 15 * mesh->numEdges, sizeof(double));
  op2_data[1] = (double *)calloc(15 * 15 * mesh->numEdges, sizeof(double));
  op_bc_data  = (double *)calloc(7 * 15 * mesh->numBoundaryEdges, sizeof(double));

  factor_data   = (double *)calloc(15 * mesh->numCells, sizeof(double));
  gFactor_data  = (double *)calloc(21 * mesh->numCells, sizeof(double));
  cFactor_data  = (double *)calloc(46 * mesh->numCells, sizeof(double));
  mmFactor_data = (double *)calloc(15 * mesh->numCells, sizeof(double));

  if(precondition) {
    glb_ind_data   = (int *)calloc(mesh->numCells, sizeof(int));
    glb_indL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
    glb_indR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
    glb_indBC_data = (int *)calloc(mesh->numBoundaryEdges, sizeof(int));
  }

  u      = op_decl_dat(mesh->cells, 15, "double", u_data, "poisson_u");
  rhs    = op_decl_dat(mesh->cells, 15, "double", rhs_data, "poisson_rhs");
  h      = op_decl_dat(mesh->cells, 1, "double", h_data, "poisson_h");
  op1    = op_decl_dat(mesh->cells, 15 * 15, "double", op1_data, "poisson_op1");
  op2[0] = op_decl_dat(mesh->edges, 15 * 15, "double", op2_data[0], "poisson_op20");
  op2[1] = op_decl_dat(mesh->edges, 15 * 15, "double", op2_data[1], "poisson_op21");
  op_bc  = op_decl_dat(mesh->bedges, 7 * 15, "double", op_bc_data, "poisson_op_bc");

  factor   = op_decl_dat(mesh->cells, 15, "double", factor_data, "poisson_factor");
  gFactor  = op_decl_dat(mesh->cells, 21, "double", gFactor_data, "poisson_gFactor");
  cFactor  = op_decl_dat(mesh->cells, 46, "double", cFactor_data, "poisson_cFactor");
  mmFactor = op_decl_dat(mesh->cells, 15, "double", mmFactor_data, "poisson_mmFactor");

  if(precondition) {
    glb_ind   = op_decl_dat(mesh->cells, 1, "int", glb_ind_data, "poisson_glb_ind");
    glb_indL  = op_decl_dat(mesh->edges, 1, "int", glb_indL_data, "poisson_glb_indL");
    glb_indR  = op_decl_dat(mesh->edges, 1, "int", glb_indR_data, "poisson_glb_indR");
    glb_indBC = op_decl_dat(mesh->bedges, 1, "int", glb_indBC_data, "poisson_glb_indBC");
  }
}

PoissonSolve::~PoissonSolve() {
  free(u_data);
  free(rhs_data);
  free(h_data);
  free(op1_data);
  free(op2_data[0]);
  free(op2_data[1]);
  free(op_bc_data);

  free(factor_data);
  free(gFactor_data);
  free(cFactor_data);
  free(mmFactor_data);

  if(precondition) {
    free(glb_ind_data);
    free(glb_indL_data);
    free(glb_indR_data);
    free(glb_indBC_data);
  }

  destroy_vec(&b);
  destroy_vec(&x);
  KSPDestroy(&ksp);
  if(matCreated)
    MatDestroy(&Amat);
}

void PoissonSolve::setDirichletBCs(int *d) {
  dirichlet[0] = d[0];
  dirichlet[1] = d[1];
  dirichlet[2] = d[2];
}

void PoissonSolve::setNeumannBCs(int *n) {
  neumann[0] = n[0];
  neumann[1] = n[1];
  neumann[2] = n[2];
}

void PoissonSolve::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double PoissonSolve::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}

void PoissonSolve::init() {
  create_vec(&b);
  create_vec(&x);
  if(precondition) {
    setGlbInd();
    op_par_loop(glb_ind_kernel, "glb_ind_kernel", mesh->edges,
                op_arg_dat(glb_ind, -2, mesh->edge2cells, 1, "int", OP_READ),
                op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));
    op_par_loop(glb_ind_kernelBC, "glb_ind_kernelBC", mesh->bedges,
                op_arg_dat(glb_ind, 0, mesh->bedge2cells, 1, "int", OP_READ),
                op_arg_dat(glb_indBC, -1, OP_ID, 1, "int", OP_WRITE));
  } else {
    create_shell_mat(&Amat);
    matCreated = true;
  }

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  KSPSetTolerances(ksp, 1e-6, 1e-50, 1e5, 1e4);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

  if(!precondition) {
    KSPSetOperators(ksp, Amat, Amat);
  }

  op_par_loop(poisson_h, "poisson_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(h, -1, OP_ID, 1, "double", OP_WRITE));
}

bool PoissonSolve::solve(op_dat b_dat, op_dat x_dat) {
  op_par_loop(poisson_apply_bc, "poisson_apply_bc", mesh->bedges,
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(b_dat, 0, mesh->bedge2cells, 15, "double", OP_INC));

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

void PoissonSolve::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  copy_u(u_d);

  op_par_loop(poisson_cells, "poisson_cells", mesh->cells,
              op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_READ),
              op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_WRITE));

  op_par_loop(poisson_edges, "poisson_edges", mesh->edges,
              op_arg_dat(u, 0, mesh->edge2cells, 15, "double", OP_READ),
              op_arg_dat(op2[0], -1, OP_ID, 15 * 15, "double", OP_READ),
              op_arg_dat(rhs, 0, mesh->edge2cells, 15, "double", OP_INC),
              op_arg_dat(u, 1, mesh->edge2cells, 15, "double", OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, 15 * 15, "double", OP_READ),
              op_arg_dat(rhs, 1, mesh->edge2cells, 15, "double", OP_INC));

  copy_rhs(rhs_d);
}

void PoissonSolve::set_op() {
  op_par_loop(poisson_op1, "poisson_op1", mesh->cells,
              op_arg_dat(mesh->cubature->J,  -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(data->Dx, -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(data->Dy, -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(cFactor, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_WRITE));

  op_par_loop(poisson_op2, "poisson_op2", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(data->mD[0], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[0], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[1], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[1], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[2], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[2], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pD[0], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pD[0], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pD[1], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pD[1], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pD[2], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pD[2], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pDy[0], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pDy[0], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pDy[1], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pDy[1], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pDy[2], 0, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->pDy[2], 1, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->edge2cells, 21, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 1, mesh->edge2cells, 21, "double", OP_READ),
              op_arg_dat(h, 0, mesh->edge2cells, 1, "double", OP_READ),
              op_arg_dat(h, 1, mesh->edge2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor, 0, mesh->edge2cells, 21, "double", OP_READ),
              op_arg_dat(gFactor, 1, mesh->edge2cells, 21, "double", OP_READ),
              op_arg_dat(factor, 0, mesh->edge2cells, 15, "double", OP_READ),
              op_arg_dat(factor, 1, mesh->edge2cells, 15, "double", OP_READ),
              op_arg_dat(op1, 0, mesh->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(op1, 1, mesh->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, 15 * 15, "double", OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, 15 * 15, "double", OP_WRITE));

  op_par_loop(poisson_op3, "poisson_op3", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(data->mD[0], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[1], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[2], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(h, 0, mesh->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(factor, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(op1, 0, mesh->bedge2cells, 15 * 15, "double", OP_INC));

  if(massMat) {
    op_par_loop(poisson_op4, "poisson_op4", mesh->cells,
                op_arg_dat(mesh->cubature->mm, -1, OP_ID, 15 * 15, "double", OP_READ),
                op_arg_dat(mmFactor, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_INC));
  }

  op_par_loop(poisson_op5, "poisson_op5", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(data->mD[0], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[1], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[2], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(h, 0, mesh->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(factor, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, 7 * 15, "double", OP_WRITE));
}

void PressureSolve::setup() {
  massMat = false;

  op_par_loop(pressure_solve_setup, "pressure_solve_setup", mesh->cells,
              op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(factor, -1, OP_ID, 15, "double", OP_WRITE));

  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(DGConstants::CUB_V), 15, factor, 0.0, cFactor);
  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), 15, factor, 0.0, gFactor);
  set_op();

  if(precondition) {
    setMatrix();
    KSPSetOperators(ksp, Amat, Amat);
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSOR);
    PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
  }
}

void ViscositySolve::setup(double mmConst) {
  massMat = true;

  op_par_loop(viscosity_solve_setup, "viscosity_solve_setup", mesh->cells,
              op_arg_dat(data->nu, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
              op_arg_gbl(&mmConst, 1, "double", OP_READ),
              op_arg_dat(factor, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(mmFactor, -1, OP_ID, 15, "double", OP_WRITE));

  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(DGConstants::CUB_V), 15, factor, 0.0, cFactor);
  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), 15, factor, 0.0, gFactor);
  set_op();

  if(precondition) {
    setMatrix();
    KSPSetOperators(ksp, Amat, Amat);
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSOR);
    PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
  }
}
