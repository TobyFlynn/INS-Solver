#ifndef __INS_2D_INS_SOLVER_BASE_H
#define __INS_2D_INS_SOLVER_BASE_H

#include "simulation_driver.h"

#include "dg_compiler_defs.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "dg_dat_pool.h"
#include "dg_linear_solvers/linear_solver.h"

#include "solvers/2d/diffusion_solver.h"

#include <string>
#include <set>

class INSSolverBase2D : public SimulationDriver {
public:
  INSSolverBase2D(DGMesh2D *m);
  INSSolverBase2D(DGMesh2D *m, const std::string &filename);
  virtual ~INSSolverBase2D() override;

  virtual void init() override;
  virtual void step() = 0;
  DG_FP get_time();
  DG_FP get_dt();

  // IO
  virtual void dump_visualisation_data(const std::string &filename) override;
  virtual void dump_checkpoint_data(const std::string &filename) override;

  // Getters (used for measurements)
  DGMesh2D* get_mesh();
  op_dat get_vel_x();
  op_dat get_vel_y();
  op_dat get_pr();

protected:
  virtual void setup_pressure_viscous_solvers(LinearSolver *pr_solver, LinearSolver *vis_solver);
  void advec_current_non_linear();
  void advec_current_non_linear_over_int();
  void advec_standard();
  void advec_standard(op_dat fx, op_dat fy, op_dat fx_old, op_dat fy_old);
  void advec_sub_cycle();
  void advec_sub_cycle(op_dat fx, op_dat fy, op_dat fx_old, op_dat fy_old);
  void advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat u_out, op_dat v_out,
                           const double t);
  void advec_sub_cycle_rhs_over_int(op_dat u_in, op_dat v_in, op_dat u_out,
                                    op_dat v_out, const double t);
  void advec_sub_cycle_rk_step(const DG_FP time_sc, const DG_FP rk_dt, op_dat u, op_dat v);
  void project_velocity_mat_mult(op_dat u, op_dat v, op_dat u_out, op_dat v_out, op_dat pen, op_dat pen_f);
  void project_velocity(op_dat dpdx, op_dat dpdy);
  DG_FP max_vel();
  void add_to_pr_history();
  void shock_capture(op_dat in0, op_dat in1);
  void filter(op_dat in);
  void zero_dat(op_dat dat);
  void update_time();
  void calc_art_vis(op_dat in, op_dat out);
  LinearSolver::Solvers set_solver_type(const std::string &str);

  DG_FP g0, a0, a1, b0, b1, dt, sub_cycle_dt, time, prev_time, h;
  int currentInd, sub_cycles;
  op_dat vel[2][2], velT[2], velTT[2], pr, n[2][2], dPdN[2], bc_types, proj_h, bc_data;

  DGMesh2D *mesh;
  bool extrapolate_initial_guess, shock_capturing, over_int_advec, gravity;
  int it_pre_sub_cycle, pr_projection_method;
  std::vector<std::pair<DG_FP,DGTempDat>> pr_history;
  std::set<std::string> values_to_save;

  // Filter params
  DG_FP filter_alpha;
  int filter_Nc, filter_sp;

  LinearSolver::Solvers pressureSolverType;
  LinearSolver::Solvers viscositySolverType;

private:
  void read_options();
  void init_dats();
  DG_FP shock_cap_calc_art_vis(op_dat in0, op_dat in1, op_dat out);

  DG_FP shock_cap_max_diff, shock_cap_smooth_tol, shock_cap_discon_tol;
  op_dat nodes_data, nodes_count, shock_cap_art_vis;
  DiffusionSolver2D *diffSolver;
};

#endif
