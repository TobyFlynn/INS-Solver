#include "measurements/2d/min_max_interface.h"

#include "op_seq.h"

#include "dg_dat_pool.h"

extern DGDatPool *dg_dat_pool;

MinMaxInterface2D::MinMaxInterface2D(INSSolverBase2D *i, const int sample_iter) : Measurement2D(i, sample_iter) {
  if(dynamic_cast<MPINSSolver2D*>(ins) == nullptr) {
    throw std::runtime_error("MinMaxInterface2D measurement can only be used with 2D multiphase solver\n");
  }

  mpins = dynamic_cast<MPINSSolver2D*>(ins);
}

void MinMaxInterface2D::measure() {
  if(!sample_this_iter())
    return;

  DGMesh2D *mesh = ins->get_mesh();
  op_dat ls = mpins->get_ls();

  DG_FP min_x = 1e6;
  DG_FP min_y = 1e6;
  DG_FP max_x = -1e6;
  DG_FP max_y = -1e6;
  op_par_loop(measure_min_max_interface, "measure_min_max_interface", mesh->cells,
              op_arg_gbl(&min_x, 1, DG_FP_STR, OP_MIN),
              op_arg_gbl(&min_y, 1, DG_FP_STR, OP_MIN),
              op_arg_gbl(&max_x, 1, DG_FP_STR, OP_MAX),
              op_arg_gbl(&max_y, 1, DG_FP_STR, OP_MAX),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(ls, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  DG_FP time = ins->get_time();
  MinMaxHistory tmp;
  tmp.time = time;
  tmp.min_x = min_x;
  tmp.min_y = min_y;
  tmp.max_x = max_x;
  tmp.max_y = max_y;
  history.push_back(tmp);
}

std::string MinMaxInterface2D::get_filename() {
  return "min_max_interface";
}

std::string MinMaxInterface2D::get_csv_header() {
  return "time,min_x,min_y,max_x,max_y";
}

std::string MinMaxInterface2D::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].min_x) + ",";
    result = result + double_to_text(history[io_count].min_y) + ",";
    result = result + double_to_text(history[io_count].max_x) + ",";
    result = result + double_to_text(history[io_count].max_y);
    io_count++;
    return result;
  } else {
    return "";
  }
}