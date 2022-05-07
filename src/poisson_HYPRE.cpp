#include "poisson_HYPRE.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

#include "dg_blas_calls.h"
#include "hypre_utils.h"

using namespace std;

extern Timing *timer;

PoissonSolveHYPRE::PoissonSolveHYPRE(DGMesh *m, INSData *nsData, LS *s) {
  mesh = m;
  data = nsData;
  ls = s;

  pMatrix = new PoissonMat(mesh);

  numberIter = 0;
  solveCount = 0;

  glb_ind_data   = (int *)calloc(mesh->numCells, sizeof(int));
  glb_indL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indBC_data = (int *)calloc(mesh->numBoundaryEdges, sizeof(int));

  orderL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  orderR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  orderBC_data = (int *)calloc(mesh->numEdges, sizeof(int));

  u_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rhs_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  in_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  out_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  pre_data = (double *)calloc(DG_NP * DG_NP * mesh->numCells, sizeof(double));

  factor_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  mmFactor_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  glb_ind   = op_decl_dat(mesh->cells, 1, "int", glb_ind_data, "poisson_glb_ind");
  glb_indL  = op_decl_dat(mesh->edges, 1, "int", glb_indL_data, "poisson_glb_indL");
  glb_indR  = op_decl_dat(mesh->edges, 1, "int", glb_indR_data, "poisson_glb_indR");
  glb_indBC = op_decl_dat(mesh->bedges, 1, "int", glb_indBC_data, "poisson_glb_indBC");

  orderL  = op_decl_dat(mesh->edges, 1, "int", orderL_data, "poisson_orderL");
  orderR  = op_decl_dat(mesh->edges, 1, "int", orderR_data, "poisson_orderR");
  orderBC = op_decl_dat(mesh->bedges, 1, "int", orderBC_data, "poisson_orderBC");

  u   = op_decl_dat(mesh->cells, DG_NP, "double", u_data, "poisson_u");
  rhs = op_decl_dat(mesh->cells, DG_NP, "double", rhs_data, "poisson_rhs");
  in  = op_decl_dat(mesh->cells, DG_NP, "double", in_data, "poisson_in");
  out = op_decl_dat(mesh->cells, DG_NP, "double", out_data, "poisson_out");
  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", pre_data, "poisson_pre");

  factor   = op_decl_dat(mesh->cells, DG_NP, "double", factor_data, "poisson_factor");
  mmFactor = op_decl_dat(mesh->cells, DG_NP, "double", mmFactor_data, "poisson_mmFactor");
}

PoissonSolveHYPRE::~PoissonSolveHYPRE() {
  free(glb_ind_data);
  free(glb_indL_data);
  free(glb_indR_data);
  free(glb_indBC_data);

  free(orderL_data);
  free(orderR_data);
  free(orderBC_data);

  free(u_data);
  free(rhs_data);
  free(in_data);
  free(out_data);
  free(pre_data);

  free(factor_data);
  free(mmFactor_data);

  HYPRE_ParCSRPCGDestroy(solver);
  HYPRE_Finalize();
}

void PoissonSolveHYPRE::init() {
  HYPREUtils::init_hypre();

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

  pMatrix->init();
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
              op_arg_dat(pMatrix->op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(b_dat,  0, mesh->bedge2cells, DG_NP, "double", OP_INC));

  HYPREUtils::dat_to_new_vec(x_dat, &ij_x, unknowns);
  HYPRE_IJVectorGetObject(ij_x, (void **) &par_x);

  HYPREUtils::dat_to_new_vec(b_dat, &ij_b, unknowns);
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

  HYPREUtils::vec_to_dat(x_dat, &ij_x, unknowns);

  HYPRE_IJVectorDestroy(ij_x);
  HYPRE_IJVectorDestroy(ij_b);
  HYPRE_IJMatrixDestroy(mat);

  return true;
}

void PoissonSolveHYPRE::setDirichletBCs(int *d) {
  pMatrix->setDirichletBCs(d);
}

void PoissonSolveHYPRE::setNeumannBCs(int *n) {
  pMatrix->setNeumannBCs(n);
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
  unknowns = mesh->get_local_vec_unknowns();
  update_glb_ind();

  op_par_loop(poisson_pr_fact, "poisson_pr_fact", mesh->cells,
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE));

  pMatrix->calc_mat(factor);
  pMatrix->transpose();

  HYPREUtils::create_matrix(&mat, unknowns);

  // Set matrix values
  setMatrix();

  HYPRE_IJMatrixAssemble(mat);
  HYPRE_IJMatrixGetObject(mat, (void **) &parcsr_mat);
}
