#include "measurements/2d/min_max_pressure.h"

#include "op_seq.h"

#include "dg_abort.h"

extern DGDatPool *dg_dat_pool;

MinMaxPressure2D::MinMaxPressure2D(INSSolverBase2D *i, const int sample_iter) : Measurement2D(i, sample_iter) {
  if(dynamic_cast<MPINSSolver2D*>(ins) == nullptr) {
    dg_abort("MinMaxInterface2D measurement can only be used with 2D multiphase solver\n");
  }

  mpins = dynamic_cast<MPINSSolver2D*>(ins);
}

void MinMaxPressure2D::measure() {
  if(!sample_this_iter())
    return;

  DGMesh2D *mesh = ins->get_mesh();
  const DG_FP alpha = mpins->get_ls_alpha();

  DG_FP min_pr = 1e18;
  DG_FP max_pr = -1e18;
  op_par_loop(measure_min_max_pressure, "measure_min_max_pressure", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&min_pr, 1, DG_FP_STR, OP_MIN),
              op_arg_gbl(&max_pr, 1, DG_FP_STR, OP_MAX),
              op_arg_dat(ins->get_pr(), -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mpins->get_ls(), -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  DG_FP time = ins->get_time();
  MinMaxHistory tmp;
  tmp.time = time;
  tmp.min_pr = min_pr;
  tmp.max_pr = max_pr;
  history.push_back(tmp);
}

std::string MinMaxPressure2D::get_filename() {
  return "min_max_pressure";
}

std::string MinMaxPressure2D::get_csv_header() {
  return "time,min_pr,max_pr,jump";
}

std::string MinMaxPressure2D::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].min_pr) + ",";
    result = result + double_to_text(history[io_count].max_pr) + ',';
    result = result + double_to_text(history[io_count].max_pr - history[io_count].min_pr);
    io_count++;
    return result;
  } else {
    return "";
  }
}