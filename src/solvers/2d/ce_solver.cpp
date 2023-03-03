#include "solvers/2d/compressible_euler_solver.h"

#include "op_seq.h"

#include <string>
#include "dg_op2_blas.h"

#include "timing.h"

extern Timing *timer;

CompressibleEulerSolver2D::CompressibleEulerSolver2D(DGMesh2D *m) {
  mesh = m;
  std::string name;
  DG_FP *tmp_np = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *tmp_g_np = (DG_FP *)calloc(DG_G_NP * mesh->cells->size, sizeof(DG_FP));
  for(int i = 0; i < 4; i++) {
    name = "Q" + std::to_string(i);
    Q[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, name.c_str());
    name = "F" + std::to_string(i);
    F[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, name.c_str());
    name = "G" + std::to_string(i);
    G[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, name.c_str());
    name = "rk_wQ" + std::to_string(i);
    rk_wQ[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, name.c_str());
    name = "rk_RHSQ0" + std::to_string(i);
    rk_RHSQ[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, name.c_str());
    name = "rk_RHSQ1" + std::to_string(i);
    rk_RHSQ[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, name.c_str());
    name = "rk_RHSQ2" + std::to_string(i);
    rk_RHSQ[2][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, name.c_str());
    name = "gQ" + std::to_string(i);
    gQ[i] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np, name.c_str());
    name = "gRHSQ" + std::to_string(i);
    gRHSQ[i] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np, name.c_str());
  }
  free(tmp_g_np);
  free(tmp_np);
}

void CompressibleEulerSolver2D::init() {
  op_par_loop(euler_2d_ic, "euler_2d_ic", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // dt = 1e-4;
  dt = 0.00125;
}

void CompressibleEulerSolver2D::step() {
  timer->startTimer("CompressibleEulerSolver2D - step");
  rhs(Q, rk_RHSQ[0]);

  op_par_loop(euler_2d_wQ_0, "euler_2d_wQ_0", mesh->cells,
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_wQ[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(rk_wQ[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(rk_wQ[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(rk_wQ[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  rhs(rk_wQ, rk_RHSQ[1]);

  op_par_loop(euler_2d_wQ_1, "euler_2d_wQ_1", mesh->cells,
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[1][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[1][3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_wQ[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(rk_wQ[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(rk_wQ[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(rk_wQ[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  rhs(rk_wQ, rk_RHSQ[2]);

  op_par_loop(euler_2d_wQ_2, "euler_2d_wQ_2", mesh->cells,
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[1][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[1][3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_RHSQ[2][3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  timer->endTimer("CompressibleEulerSolver2D - step");
}

void CompressibleEulerSolver2D::rhs(op_dat *wQ, op_dat *RHSQ) {
  op_par_loop(euler_2d_flux, "euler_2d_flux", mesh->cells,
              op_arg_dat(wQ[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(wQ[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(wQ[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(wQ[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(F[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(F[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(F[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(F[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  for(int i = 0; i < 4; i++) {
    timer->startTimer("CompressibleEulerSolver2D - div");
    mesh->div_weak(F[i], G[i], RHSQ[i]);
    timer->endTimer("CompressibleEulerSolver2D - div");
    op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, wQ[i], 0.0, gQ[i]);
  }

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(gRHSQ[0], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(gRHSQ[1], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(gRHSQ[2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(gRHSQ[3], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

  op_par_loop(euler_2d_faces, "euler_2d_faces", mesh->faces,
              op_arg_dat(mesh->order,     -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gQ[0],    -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gQ[1],    -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gQ[2],    -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gQ[3],    -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gRHSQ[0], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC),
              op_arg_dat(gRHSQ[1], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC),
              op_arg_dat(gRHSQ[2], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC),
              op_arg_dat(gRHSQ[3], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC));

  for(int i = 0; i < 4; i++) {
    op2_gemv(mesh, false, -1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, gRHSQ[i], 1.0, RHSQ[i]);
    // op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, gRHSQ[i], 1.0, RHSQ[i]);
    // inv_mass(mesh, RHSQ[i]);
  }
}

void CompressibleEulerSolver2D::dump_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(Q[0], filename.c_str());
  op_fetch_data_hdf5_file(Q[1], filename.c_str());
  op_fetch_data_hdf5_file(Q[2], filename.c_str());
  op_fetch_data_hdf5_file(Q[3], filename.c_str());
}

DG_FP CompressibleEulerSolver2D::l2_vortex_error(DG_FP time) {
  op_par_loop(euler_2d_l2_vortex_error_0, "euler_2d_l2_vortex_error_0", mesh->cells,
              op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk_wQ[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->mass(rk_wQ[0]);

  DG_FP residual = 0.0;
  op_par_loop(euler_2d_l2_vortex_error_1, "euler_2d_l2_vortex_error_1", mesh->cells,
              op_arg_gbl(&residual, 1, DG_FP_STR, OP_INC),
              op_arg_dat(rk_wQ[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  return residual;
}
