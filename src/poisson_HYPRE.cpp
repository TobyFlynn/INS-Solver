#include "poisson_HYPRE.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "dg_blas_calls.h"

using namespace std;

extern Timing *timer;

PoissonSolveHYPRE::PoissonSolveHYPRE(DGMesh *m, INSData *nsData, LS *s) {
  mesh = m;
  data = nsData;
  ls = s;

  numberIter = 0;
  solveCount = 0;

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

PoissonSolveHYPRE::~PoissonSolveHYPRE() {
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

  HYPRE_ParCSRPCGDestroy(solver);
  HYPRE_Finalize();
}

void PoissonSolveHYPRE::init() {
  HYPRE_Init();
  HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
  HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
  HYPRE_SetSpGemmUseCusparse(false);
  HYPRE_SetUseGpuRand(true);

  HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
  HYPRE_ParCSRPCGSetLogging(solver, 1);
  HYPRE_BoomerAMGCreate(&precon);
  // HYPRE_BoomerAMGSetRelaxType(precon, 3); /* 3, 4, 6, 7, 18, 11, 12 */
  // HYPRE_BoomerAMGSetRelaxOrder(precon, false); /* must be false */
  // HYPRE_BoomerAMGSetCoarsenType(precon, 8); /* 8 */
  // HYPRE_BoomerAMGSetInterpType(precon, 3); /* 3, 15, 6, 14, 18 */
  // HYPRE_BoomerAMGSetAggNumLevels(precon, 5);
  // HYPRE_BoomerAMGSetAggInterpType(precon, 5); /* 5 or 7 */
  // HYPRE_BoomerAMGSetKeepTranspose(precon, true); /* keep transpose to avoid SpMTV */
  // HYPRE_BoomerAMGSetRAP2(precon, false); /* RAP in two multiplications
  //                                           (default: FALSE) */

  op_par_loop(poisson_h, "poisson_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(h, -1, OP_ID, 1, "double", OP_WRITE));
}

void PoissonSolveHYPRE::update_glb_ind() {
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

bool PoissonSolveHYPRE::solve(op_dat b_dat, op_dat x_dat) {
  op_par_loop(poisson_apply_bc, "poisson_apply_bc", mesh->bedges,
              op_arg_dat(mesh->order,     0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(b_dat,  0, mesh->bedge2cells, DG_NP, "double", OP_INC));

  HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, unknowns - 1, &ij_x);
  HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(ij_x);

  set_x(x_dat);

  HYPRE_IJVectorAssemble(ij_x);
  HYPRE_IJVectorGetObject(ij_x, (void **) &par_x);

  HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, unknowns - 1, &ij_b);
  HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(ij_b);

  set_b(b_dat);

  HYPRE_IJVectorAssemble(ij_b);
  HYPRE_IJVectorGetObject(ij_b, (void **) &par_b);

  HYPRE_ParCSRPCGSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup,
                            precon);
  HYPRE_ParCSRPCGSetup(solver, parcsr_mat, par_b, par_x);

  timer->startKSPSolve();
  HYPRE_ParCSRPCGSolve(solver, parcsr_mat, par_b, par_x);
  timer->endKSPSolve();

  int numIt = 0;
  HYPRE_ParCSRPCGGetNumIterations(solver, &numIt);
  cout << "Number of it: " << numIt << endl;

  get_x(x_dat);

  HYPRE_IJVectorDestroy(ij_x);
  HYPRE_IJVectorDestroy(ij_b);

  return true;
}

/*
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
*/

void PoissonSolveHYPRE::set_op() {
  calc_cub_sub_mat();

  calc_gauss_sub_mat();
}

void PoissonSolveHYPRE::setDirichletBCs(int *d) {
  dirichlet[0] = d[0];
  dirichlet[1] = d[1];
  dirichlet[2] = d[2];
}

void PoissonSolveHYPRE::setNeumannBCs(int *n) {
  neumann[0] = n[0];
  neumann[1] = n[1];
  neumann[2] = n[2];
}

void PoissonSolveHYPRE::setBCValues(op_dat bc) {
  bc_dat = bc;
}

/*
double PoissonSolve::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}
*/

PressureSolveHYPRE::PressureSolveHYPRE(DGMesh *m, INSData *d, LS *s) : PoissonSolveHYPRE(m, d, s) {}

void PressureSolveHYPRE::setup() {
  unknowns = get_local_unknowns();
  update_glb_ind();

  op_par_loop(poisson_pr_fact, "poisson_pr_fact", mesh->cells,
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE));

  set_op();

  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, unknowns - 1, 0, unknowns - 1, &mat);
  HYPRE_IJMatrixSetObjectType(mat, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(mat);

  // Set matrix values
  setMatrix();

  HYPRE_IJMatrixAssemble(mat);
  HYPRE_IJMatrixGetObject(mat, (void **) &parcsr_mat);
}
