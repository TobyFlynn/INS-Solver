#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "dg_blas_calls.h"

using namespace std;

PoissonSolve::PoissonSolve(DGMesh *m, INSData *nsData) {
  mesh = m;
  data = nsData;

  numberIter = 0;
  solveCount = 0;
  pMatInit = false;
  block_jacobi_pre = false;

  glb_ind_data   = (int *)calloc(mesh->numCells, sizeof(int));
  glb_indL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indBC_data = (int *)calloc(mesh->numBoundaryEdges, sizeof(int));

  op1_data    = (double *)calloc(DG_NP * DG_NP * mesh->numCells, sizeof(double));
  op2_data[0] = (double *)calloc(DG_NP * DG_NP * mesh->numEdges, sizeof(double));
  op2_data[1] = (double *)calloc(DG_NP * DG_NP * mesh->numEdges, sizeof(double));
  op_bc_data  = (double *)calloc(DG_GF_NP * DG_NP * mesh->numBoundaryEdges, sizeof(double));

  u_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rhs_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  in_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  out_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  pre_data = (double *)calloc(DG_NP * DG_NP * mesh->numCells, sizeof(double));

  glb_ind   = op_decl_dat(mesh->cells, 1, "int", glb_ind_data, "poisson_glb_ind");
  glb_indL  = op_decl_dat(mesh->edges, 1, "int", glb_indL_data, "poisson_glb_indL");
  glb_indR  = op_decl_dat(mesh->edges, 1, "int", glb_indR_data, "poisson_glb_indR");
  glb_indBC = op_decl_dat(mesh->bedges, 1, "int", glb_indBC_data, "poisson_glb_indBC");

  op1    = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", op1_data, "poisson_op1");
  op2[0] = op_decl_dat(mesh->edges, DG_NP * DG_NP, "double", op2_data[0], "poisson_op20");
  op2[1] = op_decl_dat(mesh->edges, DG_NP * DG_NP, "double", op2_data[1], "poisson_op21");
  op_bc  = op_decl_dat(mesh->bedges, DG_GF_NP * DG_NP, "double", op_bc_data, "poisson_op_bc");

  u   = op_decl_dat(mesh->cells, DG_NP, "double", u_data, "poisson_u");
  rhs = op_decl_dat(mesh->cells, DG_NP, "double", rhs_data, "poisson_rhs");
  in  = op_decl_dat(mesh->cells, DG_NP, "double", in_data, "poisson_in");
  out = op_decl_dat(mesh->cells, DG_NP, "double", out_data, "poisson_out");
  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", pre_data, "poisson_pre");
}

PoissonSolve::~PoissonSolve() {
  free(glb_ind_data);
  free(glb_indL_data);
  free(glb_indR_data);
  free(glb_indBC_data);

  free(op1_data);
  free(op2_data[0]);
  free(op2_data[1]);
  free(op_bc_data);

  free(u_data);
  free(rhs_data);
  free(in_data);
  free(out_data);
  free(pre_data);

  if(pMatInit)
    MatDestroy(&pMat);
}

void PoissonSolve::init() {
  setGlbInd();
  op_par_loop(glb_ind_kernel, "glb_ind_kernel", mesh->edges,
              op_arg_dat(glb_ind, -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));
  op_par_loop(glb_ind_kernelBC, "glb_ind_kernelBC", mesh->bedges,
              op_arg_dat(glb_ind, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indBC, -1, OP_ID, 1, "int", OP_WRITE));
}

bool PoissonSolve::solve(op_dat b_dat, op_dat x_dat) {
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  KSPSetTolerances(ksp, 1e-6, 1e-50, 1e5, 1e4);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

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

  destroy_vec(&b);
  destroy_vec(&x);

  KSPDestroy(&ksp);

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
    op_par_loop(poisson_mm, "poisson_mm", mesh->cells,
                op_arg_gbl(&massFactor, 1, "double", OP_READ),
                op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
                op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_INC));
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

PressureSolve::PressureSolve(DGMesh *m, INSData *d) : PoissonSolve(m, d) {}

ViscositySolve::ViscositySolve(DGMesh *m, INSData *d) : PoissonSolve(m, d) {}

void PressureSolve::setup() {
  unknowns = 0;
  op_par_loop(poisson_get_num_unknowns, "poisson_get_num_unknowns", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&unknowns, 1, "int", OP_INC));

  massMat = false;
  set_op();

  // setMatrix();

  create_shell_mat(&pMat);
  pMatInit = true;
}

void ViscositySolve::setup(double mmConst) {
  unknowns = 0;
  op_par_loop(poisson_get_num_unknowns, "poisson_get_num_unknowns", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&unknowns, 1, "int", OP_INC));

  massMat = true;
  massFactor = mmConst;
  // block_jacobi_pre = true;

  set_op();

  create_shell_mat(&pMat);
  pMatInit = true;

  // PC pc;
  // KSPGetPC(ksp, &pc);
  // PCSetType(pc, PCSHELL);
  // set_shell_pc(pc);
}
