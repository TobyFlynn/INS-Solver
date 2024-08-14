#include "drivers/2d/compressible_euler_driver.h"

#include "op_seq.h"

#include "config.h"

extern Config *config;

CompressibleEulerDriver2D::CompressibleEulerDriver2D(DGMesh2D *m) {
  mesh = m;
  euler_solver = new CompressibleEuler2D(mesh);
}

CompressibleEulerDriver2D::~CompressibleEulerDriver2D() {
  delete euler_solver;
}

void CompressibleEulerDriver2D::init() {
  euler_solver->init();

  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;

  double tmp_dt = -1.0;
  config->getDouble("solver-options", "force_dt", tmp_dt);
  if(tmp_dt > 0.0)
    dt = tmp_dt;
  else
    dt = 0.5 * h / (DG_FP)((DG_ORDER + 1) * (DG_ORDER + 1));

  euler_solver->set_dt(dt);
}

void CompressibleEulerDriver2D::step() {
  euler_solver->step();
}

void CompressibleEulerDriver2D::dump_visualisation_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(euler_solver->Q[0], filename.c_str());
  op_fetch_data_hdf5_file(euler_solver->Q[1], filename.c_str());
  op_fetch_data_hdf5_file(euler_solver->Q[2], filename.c_str());
  op_fetch_data_hdf5_file(euler_solver->Q[3], filename.c_str());
}

void CompressibleEulerDriver2D::dump_checkpoint_data(const std::string &filename) {

}

DG_FP CompressibleEulerDriver2D::get_time() {
  return euler_solver->get_time();
}
  
DGMesh2D* CompressibleEulerDriver2D::get_mesh() {
  return mesh;
}

op_dat CompressibleEulerDriver2D::get_rho() {
  return euler_solver->Q[0];
}