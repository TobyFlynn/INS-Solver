#include "drivers/3d/ls_driver.h"

#include "op_seq.h"

#include "config.h"

extern Config *config;

LSDriver3D::LSDriver3D(DGMesh3D *m) {
  mesh = m;
  lsSolver = new LevelSetSolver3D(mesh);

  u = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_vel00");
  v = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_vel01");
  w = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_vel02");

  bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_bc_types");

  time = 0.0;
}

LSDriver3D::~LSDriver3D() {
  delete lsSolver;
}

void LSDriver3D::init() {
  lsSolver->init();

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

  op_par_loop(ins_3d_bc_types, "ins_3d_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  
  lsSolver->set_bc_types(bc_types);
}

void LSDriver3D::step() {
  op_par_loop(ins_3d_set_ls_vel, "ins_3d_set_ls_vel", mesh->cells,
              op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->z, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(w, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  lsSolver->step(u, v, w, dt, 1);

  time += dt;
}

void LSDriver3D::dump_visualisation_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(mesh->z, filename.c_str());
  op_fetch_data_hdf5_file(lsSolver->s, filename.c_str());
  op_fetch_data_hdf5_file(lsSolver->kink, filename.c_str());
}

void LSDriver3D::dump_checkpoint_data(const std::string &filename) {
  throw std::runtime_error("Checkpointing not supported with level-set-only");
}

DG_FP LSDriver3D::get_time() {
  return time;
}