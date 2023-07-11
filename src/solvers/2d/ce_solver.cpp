#include "solvers/2d/compressible_euler_solver.h"

#include "op_seq.h"

#include <string>
#include <stdexcept>
#include "dg_op2_blas.h"
#include "dg_dat_pool.h"

#include "timing.h"

#define L2_FREQUENCY 10

extern Timing *timer;
extern DGDatPool *dg_dat_pool;

CompressibleEulerSolver2D::CompressibleEulerSolver2D(DGMesh2D *m) {
  mesh = m;
  std::string name;
  for(int i = 0; i < 4; i++) {
    name = "Q" + std::to_string(i);
    Q[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "F" + std::to_string(i);
    F[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "G" + std::to_string(i);
    G[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "rk_wQ" + std::to_string(i);
    rk_wQ[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "rk_RHSQ0" + std::to_string(i);
    rk_RHSQ[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "rk_RHSQ1" + std::to_string(i);
    rk_RHSQ[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "rk_RHSQ2" + std::to_string(i);
    rk_RHSQ[2][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "gQ" + std::to_string(i);
  }

  l2_counter = 0;
  time = 0.0;
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

  time += dt;
  record_l2_err();
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
  }

  DGTempDat flux[4];
  flux[0] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  flux[1] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  flux[2] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  flux[3] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_2, "zero_npf_2", mesh->cells,
              op_arg_dat(flux[0].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(flux[1].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(zero_npf_2, "zero_npf_2", mesh->cells,
              op_arg_dat(flux[2].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(flux[3].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(euler_2d_faces, "euler_2d_faces", mesh->faces,
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(wQ[0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(wQ[1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(wQ[2], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(wQ[3], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(flux[0].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(flux[1].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(flux[2].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(flux[3].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  for(int i = 0; i < 4; i++) {
    op2_gemv(mesh, false, -1.0, DGConstants::LIFT, flux[i].dat, 1.0, RHSQ[i]);
  }

  dg_dat_pool->releaseTempDatCells(flux[0]);
  dg_dat_pool->releaseTempDatCells(flux[1]);
  dg_dat_pool->releaseTempDatCells(flux[2]);
  dg_dat_pool->releaseTempDatCells(flux[3]);
}

void CompressibleEulerSolver2D::dump_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(Q[0], filename.c_str());
  op_fetch_data_hdf5_file(Q[1], filename.c_str());
  op_fetch_data_hdf5_file(Q[2], filename.c_str());
  op_fetch_data_hdf5_file(Q[3], filename.c_str());
}

DG_FP CompressibleEulerSolver2D::l2_vortex_error() {
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
void CompressibleEulerSolver2D::record_l2_err() {
  if(l2_counter % L2_FREQUENCY == 0) {
    l2_err_history.push_back({time, l2_vortex_error()});
  }
  l2_counter++;
}

#include <fstream>
#include <iostream>

#include <iomanip>
#include <sstream>

std::string doubleToTextE(const double &d) {
    std::stringstream ss;
    ss << std::setprecision(15);
    ss << d;
    return ss.str();
}

#ifdef INS_MPI
#include "mpi.h"
#endif

void CompressibleEulerSolver2D::save_l2_err_history(const std::string &filename) {
  #ifdef INS_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
  #endif
  std::ofstream file(filename);

  file << "time,l2_err_pr" << std::endl;

  for(int i = 0; i < l2_err_history.size(); i++) {
    file << doubleToTextE(l2_err_history[i].first) << ",";
    file << doubleToTextE(l2_err_history[i].second) << std::endl;
  }

  file.close();
  #ifdef INS_MPI
  }
  #endif
}
