#include "drivers/2d/ls_driver.h"

#include "op_seq.h"

#include "config.h"

extern Config *config;

LSDriver2D::LSDriver2D(DGMesh2D *m) {
  mesh = m;
  lsSolver = new LevelSetSolver2D(mesh);

  u = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_vel00");
  v = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_vel01");

  bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_bc_types");

  time = 0.0;
}

LSDriver2D::~LSDriver2D() {
  delete lsSolver;
}

void LSDriver2D::init() {
  lsSolver->init();

  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  dt = 0.5 * h / (DG_FP)(DG_ORDER * DG_ORDER);

  op_par_loop(ins_bc_types, "ins_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -2, mesh->bface2nodes, 2, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  
  lsSolver->set_bc_types(bc_types);
}

void LSDriver2D::step() {
  op_par_loop(ins_2d_set_ls_vel, "ins_2d_set_ls_vel", mesh->cells,
              op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  
  lsSolver->setVelField(u, v);
  lsSolver->step(dt, 1);
}

void LSDriver2D::dump_visualisation_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(lsSolver->s, filename.c_str());

  op_dat err = op_decl_dat_temp(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "err");
  DG_FP l_max = 0.0;
  DG_FP err_dist = 5.0;
  config->getDouble("level-set-options", "err_h_dist", err_dist);
  DG_FP err_h = h;
  config->getDouble("level-set-options", "h", err_h);
  err_dist *= err_h;
  op_par_loop(ls_driver_get_err, "ls_driver_get_err", mesh->cells,
              op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&err_dist, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(lsSolver->s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(err, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_gbl(&l_max, 1, DG_FP_STR, OP_MAX));
  
  op_fetch_data_hdf5_file(err, filename.c_str());

  op_printf("L_MAX = %g\n", l_max);
  mesh->mass(err);
  DG_FP l1 = 0.0;
  op_par_loop(measure_l2_error_vortex_1, "measure_l2_error_vortex_1", mesh->cells,
              op_arg_gbl(&l1, 1, DG_FP_STR, OP_INC),
              op_arg_dat(err, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  op_printf("L_1 = %g\n", l1);
  op_free_dat_temp(err);
}

void LSDriver2D::dump_checkpoint_data(const std::string &filename) {
  throw std::runtime_error("Checkpointing not supported with level-set-only");
}

DG_FP LSDriver2D::get_time() {
  return time;
}