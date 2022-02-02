#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "dg_blas_calls.h"

using namespace std;

extern Timing *timer;

PoissonSolve::PoissonSolve(DGMesh *m, INSData *nsData, LS *s) {
  mesh = m;
  data = nsData;
  ls = s;

  numberIter = 0;
  solveCount = 0;
  pMatInit = false;
  block_jacobi_pre = false;

  glb_ind_data   = (int *)calloc(mesh->numCells, sizeof(int));
  glb_indL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indBC_data = (int *)calloc(mesh->numBoundaryEdges, sizeof(int));

  orderL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  orderR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  orderBC_data = (int *)calloc(mesh->numEdges, sizeof(int));

  op1_data    = (double *)calloc(DG_NP * DG_NP * mesh->numCells, sizeof(double));
  op2_data[0] = (double *)calloc(DG_NP * DG_NP * mesh->numEdges, sizeof(double));
  op2_data[1] = (double *)calloc(DG_NP * DG_NP * mesh->numEdges, sizeof(double));
  op_bc_data  = (double *)calloc(DG_GF_NP * DG_NP * mesh->numBoundaryEdges, sizeof(double));

  u_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rhs_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  in_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  out_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  pre_data = (double *)calloc(DG_NP * DG_NP * mesh->numCells, sizeof(double));

  factor_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  gFactor_data  = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  cFactor_data  = (double *)calloc(DG_CUB_NP * mesh->numCells, sizeof(double));
  mmFactor_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  h_data        = (double *)calloc(mesh->numCells, sizeof(double));
  gDelta_data   = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));

  glb_ind   = op_decl_dat(mesh->cells, 1, "int", glb_ind_data, "poisson_glb_ind");
  glb_indL  = op_decl_dat(mesh->edges, 1, "int", glb_indL_data, "poisson_glb_indL");
  glb_indR  = op_decl_dat(mesh->edges, 1, "int", glb_indR_data, "poisson_glb_indR");
  glb_indBC = op_decl_dat(mesh->bedges, 1, "int", glb_indBC_data, "poisson_glb_indBC");

  orderL  = op_decl_dat(mesh->edges, 1, "int", orderL_data, "poisson_orderL");
  orderR  = op_decl_dat(mesh->edges, 1, "int", orderR_data, "poisson_orderR");
  orderBC = op_decl_dat(mesh->bedges, 1, "int", orderBC_data, "poisson_orderBC");

  op1    = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", op1_data, "poisson_op1");
  op2[0] = op_decl_dat(mesh->edges, DG_NP * DG_NP, "double", op2_data[0], "poisson_op20");
  op2[1] = op_decl_dat(mesh->edges, DG_NP * DG_NP, "double", op2_data[1], "poisson_op21");
  op_bc  = op_decl_dat(mesh->bedges, DG_GF_NP * DG_NP, "double", op_bc_data, "poisson_op_bc");

  u   = op_decl_dat(mesh->cells, DG_NP, "double", u_data, "poisson_u");
  rhs = op_decl_dat(mesh->cells, DG_NP, "double", rhs_data, "poisson_rhs");
  in  = op_decl_dat(mesh->cells, DG_NP, "double", in_data, "poisson_in");
  out = op_decl_dat(mesh->cells, DG_NP, "double", out_data, "poisson_out");
  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", pre_data, "poisson_pre");

  factor   = op_decl_dat(mesh->cells, DG_NP, "double", factor_data, "poisson_factor");
  gFactor  = op_decl_dat(mesh->cells, DG_G_NP, "double", gFactor_data, "poisson_gFactor");
  cFactor  = op_decl_dat(mesh->cells, DG_CUB_NP, "double", cFactor_data, "poisson_cFactor");
  mmFactor = op_decl_dat(mesh->cells, DG_NP, "double", mmFactor_data, "poisson_mmFactor");
  h        = op_decl_dat(mesh->cells, 1, "double", h_data, "poisson_h");
  gDelta   = op_decl_dat(mesh->cells, DG_G_NP, "double", gDelta_data, "poisson_gDelta");
}

PoissonSolve::~PoissonSolve() {
  free(glb_ind_data);
  free(glb_indL_data);
  free(glb_indR_data);
  free(glb_indBC_data);

  free(orderL_data);
  free(orderR_data);
  free(orderBC_data);

  free(op1_data);
  free(op2_data[0]);
  free(op2_data[1]);
  free(op_bc_data);

  free(u_data);
  free(rhs_data);
  free(in_data);
  free(out_data);
  free(pre_data);

  free(factor_data);
  free(gFactor_data);
  free(cFactor_data);
  free(mmFactor_data);
  free(h_data);
  free(gDelta_data);

  if(pMatInit)
    MatDestroy(&pMat);

  KSPDestroy(&ksp);
}

void PoissonSolve::init() {
  op_par_loop(poisson_h, "poisson_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(h, -1, OP_ID, 1, "double", OP_WRITE));

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 1e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCGAMG);
  PCGAMGSetNSmooths(pc, 4);
  PCGAMGSetSquareGraph(pc, 1);
  PCGAMGSetNlevels(pc, 20);
  PCMGSetLevels(pc, 20, NULL);
  PCMGSetCycleType(pc, PC_MG_CYCLE_W);
  PCGAMGSetRepartition(pc, PETSC_TRUE);
  PCGAMGSetReuseInterpolation(pc, PETSC_TRUE);

  // AMGX_resources_create_simple(&rsrc);
  // AMGX_matrix_create(&matrix, rsrc, AMGX_mode_dDDI);
  // AMGX_vector_create(&rhs, rsrc, AMGX_mode_dDDI);
  // AMGX_vector_create(&soln, rsrc, AMGX_mode_dDDI);
}

void PoissonSolve::update_glb_ind() {
  setGlbInd();
  op_par_loop(copy_to_edges, "copy_to_edges", mesh->edges,
              op_arg_dat(glb_ind, -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));
  op_par_loop(copy_to_bedges, "copy_to_bedges", mesh->bedges,
              op_arg_dat(glb_ind, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indBC, -1, OP_ID, 1, "int", OP_WRITE));

  op_par_loop(copy_to_edges, "copy_to_edges", mesh->edges,
              op_arg_dat(mesh->order, -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_WRITE));
  op_par_loop(copy_to_bedges, "copy_to_bedges", mesh->bedges,
              op_arg_dat(mesh->order, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(orderBC, -1, OP_ID, 1, "int", OP_WRITE));
}

bool PoissonSolve::solve(op_dat b_dat, op_dat x_dat) {
  KSPSetOperators(ksp, pMat, pMat);

  op_par_loop(poisson_apply_bc, "poisson_apply_bc", mesh->bedges,
              op_arg_dat(mesh->order,     0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(b_dat,  0, mesh->bedge2cells, DG_NP, "double", OP_INC));

  create_vec(&b);
  create_vec(&x);

  load_vec(&b, b_dat);
  load_vec(&x, x_dat);

  timer->startKSPSolve();
  KSPSolve(ksp, b, x);
  timer->endKSPSolve();

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

  destroy_vec(&b);
  destroy_vec(&x);

  return converged;
}

// Matrix-free Mat-Vec mult function
void PoissonSolve::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  copy_vec_to_dat(u, u_d);

  op_par_loop(poisson_cells, "poisson_cells", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(u,   -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(rhs, -1, OP_ID, DG_NP, "double", OP_WRITE));

  op_par_loop(poisson_edges, "poisson_edges", mesh->edges,
              op_arg_dat(mesh->order, 0, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(u,           0, mesh->edge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[0],     -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(rhs,         0, mesh->edge2cells, DG_NP, "double", OP_INC),
              op_arg_dat(mesh->order, 1, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(u,           1, mesh->edge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[1],     -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(rhs,         1, mesh->edge2cells, DG_NP, "double", OP_INC));

  copy_dat_to_vec(rhs, rhs_d);
}

// Matrix-free block-jacobi preconditioning function
void PoissonSolve::precond(const double *in_d, double *out_d) {
  copy_vec_to_dat(in, in_d);

  op_par_loop(poisson_pre, "poisson_pre", mesh->cells,
              op_arg_dat(in,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(pre, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_WRITE));

  copy_dat_to_vec(out, out_d);
}

void PoissonSolve::set_op() {
  double tol = 1e-15;

  calc_cub_sub_mat();

  calc_gauss_sub_mat();

  if(massMat) {
    calc_mm_mat();
  }

  if(block_jacobi_pre) {
    inv_blas(mesh, op1, pre);
  }
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

PressureSolve::PressureSolve(DGMesh *m, INSData *d, LS *s) : PoissonSolve(m, d, s) {
  AMGX_SAFE_CALL(AMGX_initialize());
  AMGX_SAFE_CALL(AMGX_initialize_plugins());

  AMGX_SAFE_CALL(AMGX_config_create_from_file(&config, "/home/u1717021/Code/PhD/INS-Solver/linear_solve_config/FGMRES_CLASSICAL_AGGRESSIVE_HMIS.json"));

  AMGX_resources_create_simple(&rsrc, config);
  AMGX_matrix_create(&matrix, rsrc, AMGX_mode_dDDI);
  AMGX_vector_create(&rhs, rsrc, AMGX_mode_dDDI);
  AMGX_vector_create(&soln, rsrc, AMGX_mode_dDDI);
  AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, config);
}

PressureSolve::~PressureSolve() {
  AMGX_solver_destroy(solver);
  AMGX_vector_destroy(soln);
  AMGX_vector_destroy(rhs);
  AMGX_matrix_destroy(matrix);
  AMGX_resources_destroy(rsrc);

  AMGX_finalize_plugins();
  AMGX_finalize();
}

bool PressureSolve::solve(op_dat b_dat, op_dat x_dat) {
  op_par_loop(poisson_apply_bc, "poisson_apply_bc", mesh->bedges,
              op_arg_dat(mesh->order,     0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(b_dat,  0, mesh->bedge2cells, DG_NP, "double", OP_INC));


  uploadAMGXVec(&rhs, b_dat);
  uploadAMGXVec(&soln, x_dat);

  AMGX_solver_setup(solver, matrix);
  AMGX_solver_solve(solver, rhs, soln);

  AMGX_SOLVE_STATUS status;
  AMGX_solver_get_status(solver, &status);

  int numIt;
  AMGX_solver_get_iterations_number(solver, &numIt);
  if(status != AMGX_SOLVE_SUCCESS) {
    cout << "Number of iterations for linear solver: " << numIt << endl;
    return false;
  }
  numberIter += numIt;
  solveCount++;

  downloadAMGXVec(&soln, x_dat);

  return true;
}

ViscositySolve::ViscositySolve(DGMesh *m, INSData *d, LS *s) : PoissonSolve(m, d, s) {}

void PressureSolve::setup() {
  unknowns = get_local_unknowns();
  update_glb_ind();

  massMat = false;

  op_par_loop(poisson_pr_fact, "poisson_pr_fact", mesh->cells,
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE));

  set_op();

  // setMatrix();

  op_par_loop(poisson_transpose_cells, "poisson_transpose_cells", mesh->cells,
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_RW));

  op_par_loop(poisson_transpose_edges, "poisson_transpose_edges", mesh->edges,
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_RW),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_RW));

  timer->startBuildMat();
  // setMatrix();
  setAMGXMat();
  timer->endBuildMat();

  // create_shell_mat(&pMat);
  // pMatInit = true;
}

void ViscositySolve::setup(double mmConst) {
  unknowns = get_local_unknowns();
  update_glb_ind();

  massMat = true;
  massFactor = mmConst;
  block_jacobi_pre = true;

  op_par_loop(poisson_vis_fact, "poisson_vis_fact", mesh->cells,
              op_arg_gbl(&mmConst,   1, "double", OP_READ),
              op_arg_dat(data->mu,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(mmFactor,  -1, OP_ID, DG_NP, "double", OP_WRITE));

  set_op();

  // timer->startBuildMat();
  // setMatrix();
  // timer->endBuildMat();

  create_shell_mat();

  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCSHELL);
  set_shell_pc(pc);
}
