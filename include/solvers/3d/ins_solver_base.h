#ifndef __INS_3D_INS_SOLVER_BASE_H
#define __INS_3D_INS_SOLVER_BASE_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "dg_dat_pool.h"

class INSSolverBase3D {
public:
  INSSolverBase3D(DGMesh3D *m);
  INSSolverBase3D(DGMesh3D *m, const std::string &filename);
  ~INSSolverBase3D();

  virtual void init(const DG_FP re, const DG_FP refVel);
  virtual void step() = 0;
  DG_FP get_time();
  DG_FP get_dt();

protected:
  void advec_current_non_linear();
  void advec_standard();
  void advec_sub_cycle();
  void advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat w_in,
                           op_dat u_out, op_dat v_out, op_dat w_out,
                           const double t);
  void advec_sub_cycle_rk_step(const DG_FP time_sc, op_dat u, op_dat v, op_dat w);
  void project_velocity_mat_mult(op_dat u, op_dat v, op_dat w, op_dat u_out,
                          op_dat v_out, op_dat w_out, op_dat pen, op_dat pen_f);
  void project_velocity(op_dat dpdx, op_dat dpdy, op_dat dpdz);
  DG_FP max_vel();
  void add_to_pr_history();
  void shock_capture_filter_dat(op_dat in);

  DG_FP g0, a0, a1, b0, b1, dt, sub_cycle_dt, time, h;
  int currentInd, sub_cycles;
  op_dat vel[2][3], velT[3], velTT[3], pr, n[2][3], dPdN[2], bc_types, proj_h;

  DGMesh3D *mesh;
  bool div_div_proj, extrapolate_initial_guess, shock_cap, shock_cutoff_filter;
  int it_pre_sub_cycle;
  std::vector<std::pair<DG_FP,DGTempDat>> pr_history;

  // Filter params
  DG_FP filter_max_alpha, filter_s0, filter_k, filter_c;

private:
  void read_options();
  void init_dats();
};

#endif
