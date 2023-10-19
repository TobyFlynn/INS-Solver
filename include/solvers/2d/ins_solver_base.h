#ifndef __INS_2D_INS_SOLVER_BASE_H
#define __INS_2D_INS_SOLVER_BASE_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_2d.h"
#include "dg_dat_pool.h"

class INSSolverBase2D {
public:
  INSSolverBase2D(DGMesh2D *m);
  INSSolverBase2D(DGMesh2D *m, const std::string &filename);
  ~INSSolverBase2D();

  virtual void init(const DG_FP re, const DG_FP refVel);
  virtual void step() = 0;
  DG_FP get_time();
  DG_FP get_dt();

protected:
  void advec_current_non_linear();
  void advec_current_non_linear_over_int();
  void advec_standard();
  void advec_standard(op_dat fx, op_dat fy, op_dat fx_old, op_dat fy_old);
  void advec_sub_cycle();
  void advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat u_out, op_dat v_out,
                           const double t);
  void advec_sub_cycle_rk_step(const DG_FP time_sc, op_dat u, op_dat v);
  void project_velocity_mat_mult(op_dat u, op_dat v, op_dat u_out, op_dat v_out, op_dat pen, op_dat pen_f);
  void project_velocity(op_dat dpdx, op_dat dpdy);
  DG_FP max_vel();
  void add_to_pr_history();
  void filter(op_dat in);

  DG_FP g0, a0, a1, b0, b1, dt, sub_cycle_dt, time, h;
  int currentInd, sub_cycles;
  op_dat vel[2][2], velT[2], velTT[2], pr, n[2][2], dPdN[2], bc_types, proj_h, bc_data;

  DGMesh2D *mesh;
  bool extrapolate_initial_guess, shock_cap, over_int_advec;
  int it_pre_sub_cycle, pr_projection_method;
  std::vector<std::pair<DG_FP,DGTempDat>> pr_history;

  // Filter params
  DG_FP filter_alpha;
  int filter_Nc, filter_sp;

private:
  void read_options();
  void init_dats();
};

#endif
