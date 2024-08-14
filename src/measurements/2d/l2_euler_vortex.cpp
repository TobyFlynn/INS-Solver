#include "measurements/2d/l2_euler_vortex.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "dg_abort.h"

extern DGDatPool *dg_dat_pool;
extern DGConstants *constants;

L2EulerVortex2D::L2EulerVortex2D(SimulationDriver *d, const int sample_iter) : Measurement2D(d, sample_iter) {
  if(dynamic_cast<CompressibleEulerDriver2D*>(d) == nullptr) {
    dg_abort("L2EulerVortex2D measurement can only be used with 2D compressible Euler solver\n");
  }
  euler_driver = dynamic_cast<CompressibleEulerDriver2D*>(d);
}

void L2EulerVortex2D::measure() {
  if(!sample_this_iter())
    return;

  DGMesh2D *mesh = euler_driver->get_mesh();
  op_dat rho = euler_driver->get_rho();
  DG_FP time = euler_driver->get_time();

  DGTempDat x_cub = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat y_cub = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat rho_cub = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat err_res = dg_dat_pool->requestTempDatCells(DG_NP);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, mesh->x, 0.0, x_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, mesh->y, 0.0, y_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, rho, 0.0, rho_cub.dat);

  DG_FP *cub_w_ptr = constants->get_mat_ptr(DGConstants::CUB2D_W);
  op_par_loop(measure_error_euler_vortex, "measure_error_euler_vortex", mesh->cells,
              op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(cub_w_ptr, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(x_cub.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(y_cub.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho_cub.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, true, 1.0, DGConstants::CUB2D_INTERP, rho_cub.dat, 0.0, err_res.dat);

  dg_dat_pool->releaseTempDatCells(x_cub);
  dg_dat_pool->releaseTempDatCells(y_cub);
  dg_dat_pool->releaseTempDatCells(rho_cub);

  DG_FP residual = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&residual, 1, DG_FP_STR, OP_INC),
              op_arg_dat(err_res.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(err_res);

  L2ErrorHistory tmp;
  tmp.time = time;
  tmp.err = sqrt(residual);
  history.push_back(tmp);
}

std::string L2EulerVortex2D::get_filename() {
  return "l2_error_vortex";
}

std::string L2EulerVortex2D::get_csv_header() {
  return "time,l2_error";
}

std::string L2EulerVortex2D::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].err);
    io_count++;
    return result;
  } else {
    return "";
  }
}