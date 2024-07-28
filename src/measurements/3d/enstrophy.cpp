#include "measurements/3d/enstrophy.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "dg_abort.h"

extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

Enstrophy3D::Enstrophy3D(SimulationDriver *d, const DG_FP refMu, const DG_FP refRho,
                       const DG_FP refVel, const DG_FP refLen, const DG_FP volume,
                       const int sample_iter) : Measurement3D(d, sample_iter) {
  mu = refMu;
  rho = refRho;
  vel = refVel;
  len = refLen;
  vol = volume;
  if(dynamic_cast<INSSolverBase3D*>(d) == nullptr) {
    dg_abort("Enstrophy3D measurement can only be used with 3D single or multiphase solver\n");
  }
  ins = dynamic_cast<INSSolverBase3D*>(d);
}

void Enstrophy3D::measure() {
  if(!sample_this_iter())
    return;

  EnstrophyHistory tmp;
  tmp.time = ins->get_time() * len / vel;
  tmp.enstrophy = calc_enstrophy();
  tmp.ke = calc_ke();
  history.push_back(tmp);
}

DG_FP Enstrophy3D::calc_enstrophy() {
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
  op_par_loop(measure_enstrophy, "measure_enstrophy", mesh->cells,
              op_arg_gbl(cub_w_ptr, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 10, DG_FP_STR, OP_READ),
              op_arg_dat(cub_curl[0].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cub_curl[1].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cub_curl[2].dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, true, 1.0, DGConstants::CUB3D_INTERP, cub_curl[2].dat, 0.0, curl[2].dat);

  dg_dat_pool->releaseTempDatCells(cub_curl[0]);
  dg_dat_pool->releaseTempDatCells(cub_curl[1]);
  dg_dat_pool->releaseTempDatCells(cub_curl[2]);

  DG_FP enstrophy = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&enstrophy, 1, DG_FP_STR, OP_INC),
              op_arg_dat(curl[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(curl[0]);
  dg_dat_pool->releaseTempDatCells(curl[1]);
  dg_dat_pool->releaseTempDatCells(curl[2]);

  // Multiply by mu and divide by volume.
  // Rho is constant so just single division by refRho needed.
  return mu * enstrophy / (vol * rho * vel * vel);
}

DG_FP Enstrophy3D::calc_ke() {
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
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&kinetic_energy, 1, DG_FP_STR, OP_INC),
              op_arg_dat(ke_tmp.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(ke_tmp);

  // Multiply by 0.5 and divide by volume. Rho is constant in single phase case,
  // so cancels through in equation
  return -0.5 * kinetic_energy / (vol * vel * vel);
}

std::string Enstrophy3D::get_filename() {
  return "enstrophy";
}

std::string Enstrophy3D::get_csv_header() {
  return "time,enstrophy,kinetic_energy";
}

std::string Enstrophy3D::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].enstrophy) + ",";
    result = result + double_to_text(history[io_count].ke);
    io_count++;
    return result;
  } else {
    return "";
  }
}
