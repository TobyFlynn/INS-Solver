#include "measurements/2d/ls_advec_error.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "dg_abort.h"

extern DGDatPool *dg_dat_pool;
extern DGConstants *constants;

LSAdvecError::LSAdvecError(SimulationDriver *d, const int sample_iter) : Measurement2D(d, sample_iter) {
  if(dynamic_cast<LSDriver2D*>(d) == nullptr) {
    dg_abort("LSAdvecError measurement can only be used with 2D level-set-only solver\n");
  }

  ls_driver = dynamic_cast<LSDriver2D*>(d);
}

void LSAdvecError::measure() {
  if(!sample_this_iter())
    return;

  DGMesh2D *mesh = ls_driver->get_mesh();
  op_dat surface = ls_driver->get_surface();
  DG_FP time = ls_driver->get_time();
  DG_FP h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;

  DGTempDat l1_err_cub = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat l2_err_cub = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat x_cub = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat y_cub = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat surface_cub = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, mesh->x, 0.0, x_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, mesh->y, 0.0, y_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, surface, 0.0, surface_cub.dat);

  DG_FP err_width = 10.0 * h;
  DG_FP l_max = 0.0;
  DG_FP *cub_w_ptr = constants->get_mat_ptr(DGConstants::CUB2D_W);
  op_par_loop(measure_error_ls_advec, "measure_error_ls_advec", mesh->cells,
              op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&err_width, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(cub_w_ptr, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(x_cub.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(y_cub.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(surface_cub.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l1_err_cub.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(l2_err_cub.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_gbl(&l_max, 1, DG_FP_STR, OP_MAX));
  
  dg_dat_pool->releaseTempDatCells(x_cub);
  dg_dat_pool->releaseTempDatCells(y_cub);
  dg_dat_pool->releaseTempDatCells(surface_cub);
  DGTempDat l1_err = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat l2_err = dg_dat_pool->requestTempDatCells(DG_NP);

  op2_gemv(mesh, true, 1.0, DGConstants::CUB2D_INTERP, l1_err_cub.dat, 0.0, l1_err.dat);
  op2_gemv(mesh, true, 1.0, DGConstants::CUB2D_INTERP, l2_err_cub.dat, 0.0, l2_err.dat);

  dg_dat_pool->releaseTempDatCells(l1_err_cub);
  dg_dat_pool->releaseTempDatCells(l2_err_cub);

  DG_FP l1 = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&l1, 1, DG_FP_STR, OP_INC),
              op_arg_dat(l1_err.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  
  DG_FP l2 = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&l2, 1, DG_FP_STR, OP_INC),
              op_arg_dat(l2_err.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(l1_err);
  dg_dat_pool->releaseTempDatCells(l2_err);

  ErrorHistory tmp;
  tmp.time = time;
  tmp.l1 = l1;
  tmp.l2 = l2;
  tmp.l_max = l_max;
  history.push_back(tmp);
}

std::string LSAdvecError::get_filename() {
  return "error_ls_advec";
}

std::string LSAdvecError::get_csv_header() {
  return "time,l1_error,l2_error,l_max_error";
}

std::string LSAdvecError::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].l1) + ",";
    result = result + double_to_text(history[io_count].l2) + ",";
    result = result + double_to_text(history[io_count].l_max);
    io_count++;
    return result;
  } else {
    return "";
  }
}