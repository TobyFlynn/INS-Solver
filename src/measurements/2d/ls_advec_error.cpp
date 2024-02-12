#include "measurements/2d/ls_advec_error.h"

#include "op_seq.h"

#include "dg_dat_pool.h"
#include "dg_abort.h"

extern DGDatPool *dg_dat_pool;

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

  DGTempDat err_tmp = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(measure_l2_error_ls_advec, "measure_l2_error_ls_advec", mesh->cells,
              op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(surface, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(err_tmp.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->mass(err_tmp.dat);

  DG_FP residual = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&residual, 1, DG_FP_STR, OP_INC),
              op_arg_dat(err_tmp.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(err_tmp);

  L2ErrorHistory tmp;
  tmp.time = time;
  tmp.err = sqrt(residual);
  history.push_back(tmp);
}

std::string LSAdvecError::get_filename() {
  return "l2_error_ls_advec";
}

std::string LSAdvecError::get_csv_header() {
  return "time,l2_error";
}

std::string LSAdvecError::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].err);
    io_count++;
    return result;
  } else {
    return "";
  }
}