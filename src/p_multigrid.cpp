#include "p_multigrid.h"

#include "op_seq.h"

#include <limits>
#include <random>

#include "dg_op2_blas.h"
#include "dg_constants.h"

PMultigrid::PMultigrid(DGMesh *m, INSData *insData) {
  mesh = m;
  data = insData;
  mm = false;

  pMatrix = new PoissonMat(mesh);

  for(int i = 0; i < DG_ORDER; i++) {
    tmp_dat_data[i] = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    u_dat_data[i]   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    b_dat_data[i]   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  }

  u_rhs_data    = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rhs_rhs_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  factor_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  mmFactor_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  for(int i = 0; i < DG_ORDER; i++) {
    tmp_dat[i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_dat_data[i], "tmp");
    u_dat[i]   = op_decl_dat(mesh->cells, DG_NP, "double", u_dat_data[i], "u");
    b_dat[i]   = op_decl_dat(mesh->cells, DG_NP, "double", b_dat_data[i], "b");
  }

  u_rhs    = op_decl_dat(mesh->cells, DG_NP, "double", u_rhs_data, "u_rhs");
  rhs_rhs  = op_decl_dat(mesh->cells, DG_NP, "double", rhs_rhs_data, "rhs_rhs");
  factor   = op_decl_dat(mesh->cells, DG_NP, "double", factor_data, "pm_factor");
  mmFactor = op_decl_dat(mesh->cells, DG_NP, "double", mmFactor_data, "pm_mmFactor");
}

PMultigrid::~PMultigrid() {
  for(int i = 0; i < DG_ORDER; i++) {
    free(tmp_dat_data[i]);
    free(u_dat_data[i]);
    free(b_dat_data[i]);
  }

  free(u_rhs_data);
  free(rhs_rhs_data);
  free(factor_data);
  free(mmFactor_data);

  delete pMatrix;
}

void PMultigrid::init() {
  pMatrix->init();
}

bool PMultigrid::solve(op_dat b, op_dat x) {
  if(mm) {
    pMatrix->calc_mat_mm(factor, mmFactor);
  } else {
    pMatrix->calc_mat(factor);
  }

  op_par_loop(poisson_apply_bc, "poisson_apply_bc", mesh->bedges,
              op_arg_dat(mesh->order,     0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(pMatrix->op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(bc_dat, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(b, 0, mesh->bedge2cells, DG_NP, "double", OP_INC));

  bool converged = false;

  int p = DG_ORDER;

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_dat[p-1], -1, OP_ID, DG_NP, "double", OP_WRITE));
  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(b, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(b_dat[p-1], -1, OP_ID, DG_NP, "double", OP_WRITE));

  cycle(p);

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(u_dat[p-1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(x, -1, OP_ID, DG_NP, "double", OP_WRITE));

  // TODO convergence check
  converged = true;
  return converged;
}

void PMultigrid::cycle(int p) {
  if(mm) {
    pMatrix->calc_mat_mm(factor, mmFactor);
  } else {
    pMatrix->calc_mat(factor);
  }

  if(p == 1) {
    // u = A^-1 (F)
    sub_solve(b_dat[p-1], u_dat[p-1]);

    return;
  }

  double w = (4.0 / 3.0) * (1.0 / maxEigenValue());
  std::cout << "Factor for " << p << ": " << w << std::endl;

  // Relaxation
  // u = u + R^-1 (F - Au)
  // Probably best to switch this to some sort of convergence test
  for(int i = 0; i < 10; i++) {
    pMatrix->mult(u_dat[p-1], tmp_dat[p-1]);
    op_par_loop(p_multigrid_relaxation_jacobi, "p_multigrid_relaxation_jacobi", mesh->cells,
                op_arg_gbl(&p, 1, "int", OP_READ),
                op_arg_gbl(&w, 1, "double", OP_READ),
                op_arg_dat(tmp_dat[p-1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(b_dat[p-1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(pMatrix->op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
                op_arg_dat(u_dat[p-1], -1, OP_ID, DG_NP, "double", OP_RW));
  }


  // Restriction
  int p_new = p / 2;
  // F = I^T (F - Au)
  // u = 0
  pMatrix->mult(u_dat[p-1], tmp_dat[p-1]);
  op_par_loop(p_multigrid_restriction, "p_multigrid_restriction", mesh->cells,
              op_arg_gbl(&p, 1, "int", OP_READ),
              op_arg_dat(tmp_dat[p-1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(b_dat[p-1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(b_dat[p_new-1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(u_dat[p_new-1], -1, OP_ID, DG_NP, "double", OP_WRITE));
  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(b_dat[p_new-1]);
  // TODO Check that these should be interpolated (maybe recalculated?)
  dats_to_update.push_back(factor);
  if(mm)
    dats_to_update.push_back(mmFactor);
  mesh->update_order(p_new, dats_to_update);

  cycle(p_new);

  // Prologation
  // u = u + Iu
  std::vector<op_dat> dats_to_update2;
  dats_to_update2.push_back(u_dat[p_new-1]);
  // TODO Check that these should be interpolated (maybe recalculated?)
  dats_to_update.push_back(factor);
  if(mm)
    dats_to_update.push_back(mmFactor);
  mesh->update_order(p, dats_to_update2);
  if(mm) {
    pMatrix->calc_mat_mm(factor, mmFactor);
  } else {
    pMatrix->calc_mat(factor);
  }
  op_par_loop(p_multigrid_prolongation, "p_multigrid_prolongation", mesh->cells,
              op_arg_gbl(&p, 1, "int", OP_READ),
              op_arg_dat(u_dat[p_new-1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_dat[p-1],     -1, OP_ID, DG_NP, "double", OP_RW));


  // Relaxation
  // u = u + R^-1 (F - Au)
  // Probably best to switch this to some sort of convergence test
  for(int i = 0; i < 10; i++) {
    pMatrix->mult(u_dat[p-1], tmp_dat[p-1]);
    op_par_loop(p_multigrid_relaxation_jacobi, "p_multigrid_relaxation_jacobi", mesh->cells,
                op_arg_gbl(&p, 1, "int", OP_READ),
                op_arg_gbl(&w, 1, "double", OP_READ),
                op_arg_dat(tmp_dat[p-1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(b_dat[p-1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(pMatrix->op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
                op_arg_dat(u_dat[p-1], -1, OP_ID, DG_NP, "double", OP_RW));
  }
}

void PMultigrid::sub_solve(op_dat b_dat, op_dat u_dat) {
  unknowns = get_local_unknowns();

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 1e3);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

  Mat pMat;
  create_shell_mat(&pMat);
  KSPSetOperators(ksp, pMat, pMat);

  Vec b, x;
  create_vec(&b);
  create_vec(&x);
  load_vec(&b, b_dat);
  load_vec(&x, u_dat);

  KSPSolve(ksp, b, x);

  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  std::cout << "Number of iterations for coarsest level: " << numIt << std::endl;
  std::cout << "Converge reason: " << reason << std::endl;

  Vec solution;
  KSPGetSolution(ksp, &solution);
  store_vec(&solution, u_dat);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&pMat);
  KSPDestroy(&ksp);
}

// Matrix-free Mat-Vec mult function
void PMultigrid::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  copy_vec_to_dat(u_rhs, u_d);

  pMatrix->mult(u_rhs, rhs_rhs);

  copy_dat_to_vec(rhs_rhs, rhs_d);
}

void PMultigrid::setDirichletBCs(int *d) {
  pMatrix->setDirichletBCs(d);
}

void PMultigrid::setNeumannBCs(int *n) {
  pMatrix->setNeumannBCs(n);
}

void PMultigrid::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double PMultigrid::maxEigenValue() {
  // Get approx eigenvector using power iteration
  setRandomVector(mesh->op_tmp[0]);

  for(int i = 0; i < 10; i++) {
    pMatrix->multJacobi(mesh->op_tmp[0], mesh->op_tmp[1]);

    // Normalise vector
    double norm = 0.0;
    op_par_loop(vec_norm, "vec_norm", mesh->cells,
                op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->op_tmp[1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_gbl(&norm, 1, "double", OP_INC));
    norm = sqrt(norm);
    op_par_loop(vec_normalise, "vec_normalise", mesh->cells,
                op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->op_tmp[1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_gbl(&norm, 1, "double", OP_READ),
                op_arg_dat(mesh->op_tmp[0], -1, OP_ID, DG_NP, "double", OP_WRITE));
  }

  // Calculate eigenvalue from approx eigenvector using Rayleigh quotient
  pMatrix->multJacobi(mesh->op_tmp[0], mesh->op_tmp[1]);
  double tmp0 = 0.0;
  double tmp1 = 0.0;
  op_par_loop(rayleigh_quotient, "rayleigh_quotient", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->op_tmp[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->op_tmp[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_gbl(&tmp0, 1, "double", OP_INC),
              op_arg_gbl(&tmp1, 1, "double", OP_INC));
  return tmp0 / tmp1;
}

// TODO come up with a better way of creating a random vector in OP2
void PMultigrid::setRandomVector(op_dat vec) {
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  double randVec[DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    randVec[i] = dist(eng);
  }

  op_par_loop(set_rand_vec, "set_rand_vec", mesh->cells,
              op_arg_gbl(randVec, DG_NP, "double", OP_READ),
              op_arg_dat(vec, -1, OP_ID, DG_NP, "double", OP_WRITE));
}

PressurePMultigrid::PressurePMultigrid(DGMesh *m, INSData *d) : PMultigrid(m, d) {}

void PressurePMultigrid::setup() {
  mm = false;

  op_par_loop(poisson_pr_fact, "poisson_pr_fact", mesh->cells,
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(factor,    -1, OP_ID, DG_NP, "double", OP_WRITE));
}
