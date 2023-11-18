#include "measurements/3d/enstropy.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"

extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

Enstropy3D::Enstropy3D(INSSolverBase3D *i, const DG_FP refMu, const DG_FP refRho, 
                       const DG_FP refVel, const DG_FP volume, 
                       const int sample_iter) : Measurement3D(i, sample_iter) {
  mu = refMu;
  rho = refRho;
  vel = refVel;
  vol = volume;
}

void Enstropy3D::measure() {
  if(!sample_this_iter())
    return;

  EnstropyHistory tmp;
  tmp.time = ins->get_time();
  tmp.enstropy = calc_enstropy();
  tmp.ke = calc_ke();
  history.push_back(tmp);
}

DG_FP Enstropy3D::calc_enstropy() {
  DGMesh3D *mesh = ins->get_mesh();
  op_dat u = ins->get_vel_x();
  op_dat v = ins->get_vel_y();
  op_dat w = ins->get_vel_z();

  DGTempDat curl[3];
  curl[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  curl[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  curl[2] = dg_dat_pool->requestTempDatCells(DG_NP);

  mesh->curl(u, v, w, curl[0].dat, curl[1].dat, curl[2].dat);

  DGTempDat cub_curl[3];
  cub_curl[0] = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  cub_curl[1] = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  cub_curl[2] = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, curl[0].dat, 0.0, cub_curl[0].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, curl[1].dat, 0.0, cub_curl[1].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, curl[2].dat, 0.0, cub_curl[2].dat);

  DG_FP *cub_w_ptr = constants->get_mat_ptr(DGConstants::CUB3D_W);
  op_par_loop(measure_enstrophy_0, "measure_enstrophy_0", mesh->cells,
              op_arg_gbl(cub_w_ptr, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 10, DG_FP_STR, OP_READ),
              op_arg_dat(cub_curl[0].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cub_curl[1].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cub_curl[2].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, true, 1.0, DGConstants::CUB3D_INTERP, cub_curl[2].dat, 0.0, curl[2].dat);

  dg_dat_pool->releaseTempDatCells(cub_curl[0]);
  dg_dat_pool->releaseTempDatCells(cub_curl[1]);
  dg_dat_pool->releaseTempDatCells(cub_curl[2]);

  DG_FP enstropy = 0.0;
  op_par_loop(measure_enstrophy_1, "measure_enstrophy_1", mesh->cells,
              op_arg_gbl(&enstropy, 1, DG_FP_STR, OP_INC),
              op_arg_dat(curl[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(curl[0]);
  dg_dat_pool->releaseTempDatCells(curl[1]);
  dg_dat_pool->releaseTempDatCells(curl[2]);

  // Multiply by mu and divide by volume. 
  // Rho is constant so just single division by refRho needed.
  return mu * enstropy / (vol * rho * vel * vel);
}

DG_FP Enstropy3D::calc_ke() {
  DGMesh3D *mesh = ins->get_mesh();
  op_dat u = ins->get_vel_x();
  op_dat v = ins->get_vel_y();
  op_dat w = ins->get_vel_z();

  DGTempDat cub_ke_tmp = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);

  DGTempDat cub_vel[3];
  cub_vel[0] = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  cub_vel[1] = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  cub_vel[2] = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, u, 0.0, cub_vel[0].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, v, 0.0, cub_vel[1].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, w, 0.0, cub_vel[2].dat);

  DG_FP *cub_w_ptr = constants->get_mat_ptr(DGConstants::CUB3D_W);
  op_par_loop(measure_ke, "measure_ke", mesh->cells,
              op_arg_gbl(cub_w_ptr, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 10, DG_FP_STR, OP_READ),
              op_arg_dat(cub_vel[0].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cub_vel[1].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cub_vel[2].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cub_ke_tmp.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_WRITE));

  dg_dat_pool->releaseTempDatCells(cub_vel[0]);
  dg_dat_pool->releaseTempDatCells(cub_vel[1]);
  dg_dat_pool->releaseTempDatCells(cub_vel[2]);

  DGTempDat ke_tmp = dg_dat_pool->requestTempDatCells(DG_NP);

  op2_gemv(mesh, true, 1.0, DGConstants::CUB3D_INTERP, cub_ke_tmp.dat, 0.0, ke_tmp.dat);

  dg_dat_pool->releaseTempDatCells(cub_ke_tmp);

  DG_FP kinetic_energy = 0.0;
  op_par_loop(measure_enstrophy_1, "measure_enstrophy_1", mesh->cells,
              op_arg_gbl(&kinetic_energy, 1, DG_FP_STR, OP_INC),
              op_arg_dat(ke_tmp.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(ke_tmp);

  // Multiply by 0.5 and divide by volume. Rho is constant in single phase case,
  // so cancels through in equation
  return 0.5 * kinetic_energy / (vol * vel * vel);
}

std::string Enstropy3D::get_filename() {
  return "enstropy";
}

std::string Enstropy3D::get_csv_header() {
  return "time,enstropy,kinetic_energy";
}

std::string Enstropy3D::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].enstropy) + ",";
    result = result + double_to_text(history[io_count].ke);
    io_count++;
    return result;
  } else {
    return "";
  }
}