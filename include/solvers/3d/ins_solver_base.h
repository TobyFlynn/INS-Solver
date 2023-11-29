#ifndef __INS_3D_INS_SOLVER_BASE_H
#define __INS_3D_INS_SOLVER_BASE_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "dg_dat_pool.h"

#include <string>
#include <set>

class INSSolverBase3D {
public:
  INSSolverBase3D(DGMesh3D *m);
  INSSolverBase3D(DGMesh3D *m, const std::string &filename);
  ~INSSolverBase3D();

  virtual void init(const DG_FP re, const DG_FP refVel);
  virtual void step() = 0;
  DG_FP get_time();
  DG_FP get_dt();

  // IO
  virtual void dump_visualisation_data(const std::string &filename);
  virtual void dump_checkpoint_data(const std::string &filename);

  // Getters (used for measurements)
  DGMesh3D* get_mesh();
  op_dat get_vel_x();
  op_dat get_vel_y();
  op_dat get_vel_z();
  op_dat get_pr();

protected:
  void advec_current_non_linear();
  void advec_current_non_linear_over_int();
  void advec_standard();
  void advec_sub_cycle();
  void advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat w_in,
                           op_dat u_out, op_dat v_out, op_dat w_out,
                           const double t);
  void advec_sub_cycle_rhs_over_int(op_dat u_in, op_dat v_in, op_dat w_in,
                           op_dat u_out, op_dat v_out, op_dat w_out,
                           const double t);
  void advec_sub_cycle_rk_step(const DG_FP time_sc, const DG_FP rk_dt, op_dat u, op_dat v, op_dat w);
  void project_velocity_mat_mult(op_dat u, op_dat v, op_dat w, op_dat u_out,
                          op_dat v_out, op_dat w_out, op_dat pen, op_dat pen_f);
  void project_velocity(op_dat dpdx, op_dat dpdy, op_dat dpdz);
  DG_FP max_vel();
  void add_to_pr_history();
  void filter(op_dat in);
  void zero_dat(op_dat dat);
  void update_time();

  DG_FP g0, a0, a1, b0, b1, dt, sub_cycle_dt, time, prev_time, h;
  int currentInd, sub_cycles;
  op_dat vel[2][3], velT[3], velTT[3], pr, n[2][3], dPdN[2], bc_types, proj_h;

  DGMesh3D *mesh;
  bool extrapolate_initial_guess;
  bool over_int_advec, filter_advec;
  int it_pre_sub_cycle, pr_projection_method;
  std::vector<std::pair<DG_FP,DGTempDat>> pr_history;
  std::set<std::string> values_to_save;

  // Filter params
  DG_FP filter_alpha;
  int filter_Nc, filter_sp;

private:
  void read_options();
  void init_dats();
};

#endif
