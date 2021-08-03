#include "poisson.h"

#include "op_seq.h"

#include <iostream>

#include "dg_constants.h"
#include "dg_blas_calls.h"
#include "dg_compiler_defs.h"

extern DGConstants *constants;

using namespace std;

PoissonSolve::PoissonSolve(DGMesh *m, INSData *d, bool p) {
  mesh = m;
  data = d;
  precondition = p;

  h_data       = (double *)calloc(mesh->numCells, sizeof(double));
  op1_data     = (double *)calloc(10 * 10 * mesh->numCells, sizeof(double));
  op2_data[0]  = (double *)calloc(10 * 10 * mesh->numEdges, sizeof(double));
  op2_data[1]  = (double *)calloc(10 * 10 * mesh->numEdges, sizeof(double));
  op_bc_data   = (double *)calloc(6 * 10 * mesh->numBoundaryEdges, sizeof(double));
  cFactor_data = (double *)calloc(36 * mesh->numCells, sizeof(double));

  if(precondition) {
    glb_ind_data   = (int *)calloc(mesh->numCells, sizeof(int));
    glb_indL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
    glb_indR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
    glb_indBC_data = (int *)calloc(mesh->numBoundaryEdges, sizeof(int));
  }

  tmp_data = (double *)calloc(10 * 10 * mesh->numCells, sizeof(double));
  pre_data = (double *)calloc(10 * 10 * mesh->numCells, sizeof(double));

  h       = op_decl_dat(mesh->cells, 1, "double", h_data, "poisson_h");
  op1     = op_decl_dat(mesh->cells, 10 * 10, "double", op1_data, "poisson_op1");
  op2[0]  = op_decl_dat(mesh->edges, 10 * 10, "double", op2_data[0], "poisson_op20");
  op2[1]  = op_decl_dat(mesh->edges, 10 * 10, "double", op2_data[1], "poisson_op21");
  op_bc   = op_decl_dat(mesh->bedges, 6 * 10, "double", op_bc_data, "poisson_op_bc");
  cFactor = op_decl_dat(mesh->cells, 36, "double", cFactor_data, "poisson_cFactor");

  if(precondition) {
    glb_ind   = op_decl_dat(mesh->cells, 1, "int", glb_ind_data, "poisson_glb_ind");
    glb_indL  = op_decl_dat(mesh->edges, 1, "int", glb_indL_data, "poisson_glb_indL");
    glb_indR  = op_decl_dat(mesh->edges, 1, "int", glb_indR_data, "poisson_glb_indR");
    glb_indBC = op_decl_dat(mesh->bedges, 1, "int", glb_indBC_data, "poisson_glb_indBC");
  }

  tmp = op_decl_dat(mesh->cells, 10 * 10, "double", tmp_data, "poisson_tmp");
  pre = op_decl_dat(mesh->cells, 10 * 10, "double", pre_data, "poisson_pre");
}

PoissonSolve::~PoissonSolve() {
  free(h_data);
  free(op1_data);
  free(op2_data[0]);
  free(op2_data[1]);
  free(op_bc_data);
  free(cFactor_data);

  if(precondition) {
    free(glb_ind_data);
    free(glb_indL_data);
    free(glb_indR_data);
    free(glb_indBC_data);
  }

  free(tmp_data);
  free(pre_data);

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
  u   = data->tmp_dg_np[4];
  rhs = data->tmp_dg_np[5];
  in  = data->tmp_dg_np[6];
  out = data->tmp_dg_np[7];

  factor   = data->tmp_dg_np[8];
  mmFactor = data->tmp_dg_np[9];
  gFactor  = data->tmp_dg_g_np[3];

  create_vec(&b, 10);
  create_vec(&x, 10);
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
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, 18, "double", OP_READ),
              op_arg_dat(b_dat,  0, mesh->bedge2cells, 10, "double", OP_INC));

  load_vec(&b, b_dat, 10);
  load_vec(&x, x_dat, 10);

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
  copy_vec_to_dat(u, u_d);

  op_par_loop(poisson_cells, "poisson_cells", mesh->cells,
              op_arg_dat(u,   -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, 10 * 10, "double", OP_READ),
              op_arg_dat(rhs, -1, OP_ID, 10, "double", OP_WRITE));

  op_par_loop(poisson_edges, "poisson_edges", mesh->edges,
              op_arg_dat(u,       0, mesh->edge2cells, 10, "double", OP_READ),
              op_arg_dat(op2[0], -1, OP_ID, 10 * 10, "double", OP_READ),
              op_arg_dat(rhs,     0, mesh->edge2cells, 10, "double", OP_INC),
              op_arg_dat(u,       1, mesh->edge2cells, 10, "double", OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, 10 * 10, "double", OP_READ),
              op_arg_dat(rhs,     1, mesh->edge2cells, 10, "double", OP_INC));

  copy_dat_to_vec(rhs, rhs_d);
}

void PoissonSolve::precond(const double *in_d, double *out_d) {
  copy_vec_to_dat(in, in_d);

  op_par_loop(poisson_pre, "poisson_pre", mesh->cells,
              op_arg_dat(in,  -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(pre, -1, OP_ID, 10 * 10, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, 10, "double", OP_WRITE));

  copy_dat_to_vec(out, out_d);
}

void PoissonSolve::set_op() {
  op_par_loop(poisson_op1, "poisson_op1", mesh->cells,
              op_arg_dat(mesh->cubature->J, -1, OP_ID, 36, "double", OP_READ),
              op_arg_dat(data->Dx, -1, OP_ID, 36 * 10, "double", OP_READ),
              op_arg_dat(data->Dy, -1, OP_ID, 36 * 10, "double", OP_READ),
              op_arg_dat(cFactor,  -1, OP_ID, 36, "double", OP_READ),
              op_arg_dat(op1,      -1, OP_ID, 10 * 10, "double", OP_WRITE));

  op_par_loop(poisson_op2, "poisson_op2", mesh->edges,
              op_arg_dat(mesh->edgeNum,  -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,  -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(data->mD[0],     0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[0],     1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[1],     0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[1],     1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[2],     0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[2],     1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pD[0],     0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pD[0],     1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pD[1],     0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pD[1],     1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pD[2],     0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pD[2],     1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pDy[0],    0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pDy[0],    1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pDy[1],    0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pDy[1],    1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pDy[2],    0, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->pDy[2],    1, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->edge2cells, 18, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 1, mesh->edge2cells, 18, "double", OP_READ),
              op_arg_dat(h, 0, mesh->edge2cells, 1, "double", OP_READ),
              op_arg_dat(h, 1, mesh->edge2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor, 0, mesh->edge2cells, 18, "double", OP_READ),
              op_arg_dat(gFactor, 1, mesh->edge2cells, 18, "double", OP_READ),
              op_arg_dat(factor,  0, mesh->edge2cells, 10, "double", OP_READ),
              op_arg_dat(factor,  1, mesh->edge2cells, 10, "double", OP_READ),
              op_arg_dat(op1,     0, mesh->edge2cells, 10 * 10, "double", OP_INC),
              op_arg_dat(op1,     1, mesh->edge2cells, 10 * 10, "double", OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, 10 * 10, "double", OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, 10 * 10, "double", OP_WRITE));

  op_par_loop(poisson_op3, "poisson_op3", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(data->mD[0],     0, mesh->bedge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[1],     0, mesh->bedge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[2],     0, mesh->bedge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, 18, "double", OP_READ),
              op_arg_dat(h,               0, mesh->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor,         0, mesh->bedge2cells, 18, "double", OP_READ),
              op_arg_dat(factor,          0, mesh->bedge2cells, 10, "double", OP_READ),
              op_arg_dat(op1,             0, mesh->bedge2cells, 10 * 10, "double", OP_INC));

  if(massMat) {
    op_par_loop(poisson_op4, "poisson_op4", mesh->cells,
                op_arg_dat(mesh->cubature->mm, -1, OP_ID, 10 * 10, "double", OP_READ),
                op_arg_dat(mmFactor, -1, OP_ID, 10, "double", OP_READ),
                op_arg_dat(op1,      -1, OP_ID, 10 * 10, "double", OP_INC),
                op_arg_dat(tmp,      -1, OP_ID, 10 * 10, "double", OP_WRITE));

    inv_blas(mesh, op1, pre);
  }

  op_par_loop(poisson_op5, "poisson_op5", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(data->mD[0],     0, mesh->bedge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[1],     0, mesh->bedge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(data->mD[2],     0, mesh->bedge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, 18, "double", OP_READ),
              op_arg_dat(h,               0, mesh->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor,         0, mesh->bedge2cells, 18, "double", OP_READ),
              op_arg_dat(factor,          0, mesh->bedge2cells, 10, "double", OP_READ),
              op_arg_dat(op_bc,          -1, OP_ID, 6 * 10, "double", OP_WRITE));
}

void PressureSolve::setup() {
  massMat = false;

  op_par_loop(pressure_solve_setup, "pressure_solve_setup", mesh->cells,
              op_arg_dat(data->rho, -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, 10, "double", OP_WRITE));

  op2_gemv(true, 36, 10, 1.0, constants->get_ptr(DGConstants::CUB_V), 10, factor, 0.0, cFactor);
  op2_gemv(true, 18, 10, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), 10, factor, 0.0, gFactor);
  set_op();

  if(precondition) {
    setMatrix();
    KSPSetOperators(ksp, Amat, Amat);

    PC pc;
    KSPGetPC(ksp, &pc);
    // PCSetType(pc, PCILU);

    // PCSetType(pc, PCASM);
    // PCASMSetTotalSubdomains(pc, 2, NULL, NULL);

    PCSetType(pc, PCGAMG);
    PCGAMGSetNSmooths(pc, 3);
    PCGAMGSetSquareGraph(pc, 0);
    PCGAMGSetNlevels(pc, 5);
    PCMGSetLevels(pc, 10, NULL);
    PCMGSetCycleType(pc, PC_MG_CYCLE_W);
    PCGAMGSetRepartition(pc, PETSC_TRUE);
    PCGAMGSetReuseInterpolation(pc, PETSC_TRUE);

    // PetscOptionsSetValue(NULL, "-pc_type", "gamg");
    // PetscOptionsSetValue(NULL, "-pc_gamg_agg_nsmooths", "3");
    // PetscOptionsSetValue(NULL, "-pc_mg_cycle_type", "v");
    // PetscOptionsSetValue(NULL, "-mg_levels_ksp_max_it", "10");

    // PetscOptionsSetValue(NULL, "-pc_type", "ml");
    // PetscOptionsSetValue(NULL, "-pc_ml_maxNlevels", "10");
    // PetscOptionsSetValue(NULL, "-pc_mg_cycle_type", "w");
    // PetscOptionsSetValue(NULL, "-mg_levels_ksp_max_it", "4");

    // PetscOptionsSetValue(NULL, "-pc_type", "hypre");
    // PetscOptionsSetValue(NULL, "-pc_hypre_type", "boomeramg");
    // PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_coarsen_type", "PMIS");
    // PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_interp_type", "direct");
    // PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type", "W");
    // PetscOptionsSetValue(NULL, "--pc_hypre_boomeramg_agg_nl", "1");

    // PetscOptionsSetValue(NULL, "-help", "");
    // PCSetFromOptions(pc);

    // PCSetType(pc, PCML);

    // PCSetType(pc, PCHYPRE);

    KSPSetOperators(ksp, Amat, Amat);
  }
}

void ViscositySolve::setup(double mmConst) {
  massMat = true;

  if(!matCreated) {
    create_shell_mat(&Amat);
    matCreated = true;
    KSPSetOperators(ksp, Amat, Amat);
  }

  op_par_loop(viscosity_solve_setup, "viscosity_solve_setup", mesh->cells,
              op_arg_dat(data->nu,  -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, 10, "double", OP_READ),
              op_arg_gbl(&mmConst,   1, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, 10, "double", OP_WRITE),
              op_arg_dat(mmFactor,  -1, OP_ID, 10, "double", OP_WRITE));

  op2_gemv(true, 36, 10, 1.0, constants->get_ptr(DGConstants::CUB_V), 10, factor, 0.0, cFactor);
  op2_gemv(true, 18, 10, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), 10, factor, 0.0, gFactor);
  set_op();

  if(precondition) {
    setMatrix();
    KSPSetOperators(ksp, Amat, Amat);

    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSHELL);
    set_shell_pc(pc);
  }
}
