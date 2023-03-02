#include "linear_solvers/pmultigrid.h"

#include "op_seq.h"

#include <random>
#include <string>

#include "matrices/3d/poisson_semi_matrix_free_3d.h"
#include "utils.h"
#include "timing.h"

extern Timing *timer;

#define RAND_VEC_SIZE 25

PMultigridPoissonSolver::PMultigridPoissonSolver(DGMesh *m) {
  mesh = m;
  DG_FP *tmp_data = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  std::string name;
  for(int i = 0; i < DG_ORDER; i++) {
    name = "p_multigrid_tmp" + std::to_string(i);
    tmp_dat[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, name.c_str());
    name = "p_multigrid_u" + std::to_string(i);
    u_dat[i]   = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, name.c_str());
    name = "p_multigrid_b" + std::to_string(i);
    b_dat[i]   = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, name.c_str());
  }
  eg_tmp_0 = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, "p_multigrid_eg_tmp_0");
  eg_tmp_1 = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, "p_multigrid_eg_tmp_1");

  rk[0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, "p_multigrid_rk0");
  rk[1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, "p_multigrid_rk1");
  rk[2] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, "p_multigrid_rk2");
  rkQ   = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_data, "p_multigrid_rkQ");
  free(tmp_data);

  coarseSolver = new PETScAMGSolver(mesh);
}

PMultigridPoissonSolver::~PMultigridPoissonSolver() {
  delete coarseSolver;
}

bool PMultigridPoissonSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("PMultigridPoissonSolver - solve");
  if(bc)
    matrix->apply_bc(rhs, bc);

  int order = DG_ORDER;

  setupDirectSolve();

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(ans, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u_dat[order - 1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(rhs, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(b_dat[order - 1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  cycle(order);

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(u_dat[order - 1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(ans, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  timer->endTimer("PMultigridPoissonSolver - solve");

  return true;
}

void PMultigridPoissonSolver::cycle(int order) {
  if(order == 1) {
    // u = A^-1 (F)
    if(coarseMatCalcRequired) {
      timer->startTimer("PMultigridPoissonSolver - Calc Mat");
      coarseMatrix->calc_mat();
      coarseMatCalcRequired = false;
      timer->endTimer("PMultigridPoissonSolver - Calc Mat");
    }

    timer->startTimer("PMultigridPoissonSolver - Direct Solve");
    coarseSolver->solve(b_dat[order-1], u_dat[order-1]);
    timer->endTimer("PMultigridPoissonSolver - Direct Solve");
    return;
  }
  PoissonSemiMatrixFree3D *smf_mat = (PoissonSemiMatrixFree3D *)matrix;
  timer->startTimer("PMultigridPoissonSolver - Calc Mat Partial");
  smf_mat->calc_mat_partial();
  timer->endTimer("PMultigridPoissonSolver - Calc Mat Partial");

  // Calc factor for relaxation
  timer->startTimer("PMultigridPoissonSolver - w");
  w = (4.0 / 3.0) * (1.0 / maxEigenValue());
  timer->endTimer("PMultigridPoissonSolver - w");

  // Relaxation
  // u = u + R^-1 (F - Au)
  timer->startTimer("PMultigridPoissonSolver - Relaxation");
  for(int i = 0; i < num_pre_relax_iter; i++) {
    smoother(order);
  }
  timer->endTimer("PMultigridPoissonSolver - Relaxation");

  // Restriction
  timer->startTimer("PMultigridPoissonSolver - Restriction");
  int order_new = order / 2;
  // F = I^T (F - Au)
  matrix->mult(u_dat[order-1], tmp_dat[order-1]);
  op_par_loop(p_multigrid_restriction, "p_multigrid_restriction", mesh->cells,
              op_arg_dat(mesh->order,        -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(tmp_dat[order-1],   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(b_dat[order-1],     -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(b_dat[order_new-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(u_dat[order_new-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(b_dat[order_new-1]);
  timer->startTimer("PMultigridPoissonSolver - Interp");
  mesh->update_order(order_new, dats_to_update);
  timer->endTimer("PMultigridPoissonSolver - Interp");
  timer->endTimer("PMultigridPoissonSolver - Restriction");

  cycle(order_new);

  // Prologation
  // u = u + Iu
  timer->startTimer("PMultigridPoissonSolver - Prolongation");
  std::vector<op_dat> dats_to_update2;
  dats_to_update2.push_back(u_dat[order_new-1]);
  timer->startTimer("PMultigridPoissonSolver - Interp");
  mesh->update_order(order, dats_to_update2);
  timer->endTimer("PMultigridPoissonSolver - Interp");

  op_par_loop(p_multigrid_prolongation, "p_multigrid_prolongation", mesh->cells,
              op_arg_dat(mesh->order,        -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(u_dat[order_new-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u_dat[order-1],     -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("PMultigridPoissonSolver - Prolongation");

  timer->startTimer("PMultigridPoissonSolver - Calc Mat Partial");
  smf_mat->calc_mat_partial();
  timer->endTimer("PMultigridPoissonSolver - Calc Mat Partial");

  // Relaxation
  // u = u + R^-1 (F - Au)
  timer->startTimer("PMultigridPoissonSolver - Relaxation");
  for(int i = 0; i < num_post_relax_iter; i++) {
    smoother(order);
  }
  timer->endTimer("PMultigridPoissonSolver - Relaxation");
}

void PMultigridPoissonSolver::set_coarse_matrix(PoissonCoarseMatrix *c_mat) {
  coarseMatrix = c_mat;
  coarseSolver->set_matrix(coarseMatrix);
  coarseMatCalcRequired = true;
}

void PMultigridPoissonSolver::setupDirectSolve() {
  coarseSolver->set_bcs(bc);
  coarseSolver->set_nullspace(nullspace);
  coarseSolver->set_tol(1e-2);
}

DG_FP PMultigridPoissonSolver::maxEigenValue() {
  // Get approx eigenvector using power iteration
  timer->startTimer("PMultigridPoissonSolver - Random Vec");
  setRandomVector(eg_tmp_0);
  timer->endTimer("PMultigridPoissonSolver - Random Vec");

  for(int i = 0; i < 10; i++) {
    matrix->multJacobi(eg_tmp_0, eg_tmp_1);

    // Normalise vector
    DG_FP norm = 0.0;
    op_par_loop(p_multigrid_vec_norm, "p_multigrid_vec_norm", mesh->cells,
                op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(eg_tmp_1, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(&norm, 1, DG_FP_STR, OP_INC));

    norm = sqrt(norm);
    op_par_loop(p_multigrid_vec_normalise, "p_multigrid_vec_normalise", mesh->cells,
                op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(eg_tmp_1, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(&norm, 1, DG_FP_STR, OP_READ),
                op_arg_dat(eg_tmp_0, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  // Calculate eigenvalue from approx eigenvector using Rayleigh quotient
  matrix->multJacobi(eg_tmp_0, eg_tmp_1);
  DG_FP tmp0 = 0.0;
  DG_FP tmp1 = 0.0;
  op_par_loop(p_multigrid_rayleigh_quotient, "p_multigrid_rayleigh_quotient", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(eg_tmp_0, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(eg_tmp_1, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(&tmp0, 1, DG_FP_STR, OP_INC),
              op_arg_gbl(&tmp1, 1, DG_FP_STR, OP_INC));
  return tmp0 / tmp1;
}

void PMultigridPoissonSolver::setRandomVector(op_dat vec) {
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<DG_FP> dist(0.0, 1.0);

  DG_FP rand_vec[RAND_VEC_SIZE];
  for(int i = 0; i < RAND_VEC_SIZE; i++) {
    rand_vec[i] = dist(eng);
  }

  DG_FP *vec_ptr = getOP2PtrHost(vec, OP_WRITE);

  #pragma omp parallel for
  for(int i = 0; i < vec->set->size * vec->dim; i++) {
    vec_ptr[i] = rand_vec[i % RAND_VEC_SIZE];
  }

  releaseOP2PtrHost(vec, OP_WRITE, vec_ptr);
}

void PMultigridPoissonSolver::smoother(const int order) {
  matrix->mult(u_dat[order-1], tmp_dat[order-1]);
  op_par_loop(p_multigrid_relaxation_jacobi, "p_multigrid_relaxation_jacobi", mesh->cells,
              op_arg_gbl(&w, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->order,      -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(tmp_dat[order-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(b_dat[order-1],   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(matrix->op1,      -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u_dat[order-1],   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
}
/*
void PMultigridPoissonSolver::smoother(const int order) {
  int x = -1;
  const double dt = 1e-4;
  op_par_loop(p_multigrid_runge_kutta_0, "p_multigrid_runge_kutta_0", mesh->cells,
              op_arg_gbl(&x,     1, "int", OP_READ),
              op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
              op_arg_dat(u_dat[order-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(b_dat[order-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rkQ,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  for(int j = 0; j < 3; j++) {
    matrix->mult(rkQ, rk[j]);

    if(j != 2) {
      op_par_loop(p_multigrid_runge_kutta_0, "p_multigrid_runge_kutta_0", mesh->cells,
                  op_arg_gbl(&j,     1, "int", OP_READ),
                  op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
                  op_arg_dat(u_dat[order-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(b_dat[order-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rk[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rk[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rkQ,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
    }
  }

  op_par_loop(p_multigrid_runge_kutta_1, "p_multigrid_runge_kutta_1", mesh->cells,
              op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
              op_arg_dat(u_dat[order-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(rk[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(b_dat[order-1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
}
*/
