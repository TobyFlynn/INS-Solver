#include "measurements/2d/l2_error_vortex.h"

#include "op_seq.h"

#include <fstream>
#include <iostream>

#include "dg_dat_pool.h"

extern DGDatPool *dg_dat_pool;

L2ErrorVortex2D::L2ErrorVortex2D(INSSolverBase2D *i, const int sample_iter) : Measurement2D(i, sample_iter) {
  
}

void L2ErrorVortex2D::measure() {
  if(!sample_this_iter())
    return;

  DGMesh2D *mesh = ins->get_mesh();
  op_dat u = ins->get_vel_x();
  DG_FP time = ins->get_time();

  DGTempDat err_tmp = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(measure_l2_error_vortex_0, "ins_2d_l2_vortex_error_0", mesh->cells,
              op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(err_tmp.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->mass(err_tmp.dat);

  DG_FP residual = 0.0;
  op_par_loop(measure_l2_error_vortex_1, "measure_l2_error_vortex_1", mesh->cells,
              op_arg_gbl(&residual, 1, DG_FP_STR, OP_INC),
              op_arg_dat(err_tmp.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(err_tmp);

  L2ErrorHistory tmp;
  tmp.time = time;
  tmp.err = sqrt(residual);
  history.push_back(tmp);
}

void L2ErrorVortex2D::output(std::string &path) {
  std::ofstream file(path + "l2_error_vortex.txt");

  file << "time,l2_error" << std::endl;
  for(auto &sample : history) {
    file << double_to_text(sample.time) << "," << double_to_text(sample.err) << std::endl;
  }

  file.close();
}