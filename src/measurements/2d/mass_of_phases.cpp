#include "measurements/2d/mass_of_phases.h"

#include "op_seq.h"

#include "dg_dat_pool.h"
#include "dg_abort.h"

extern DGDatPool *dg_dat_pool;

MassOfPhases2D::MassOfPhases2D(SimulationDriver *d, const int sample_iter) : Measurement2D(d, sample_iter) {
  if(dynamic_cast<MPINSSolver2D*>(d) == nullptr) {
    dg_abort("MassOfPhases2D measurement can only be used with 2D multiphase solver\n");
  }

  mpins = dynamic_cast<MPINSSolver2D*>(d);
}

void MassOfPhases2D::measure() {
  if(!sample_this_iter())
    return;

  DGMesh2D *mesh = mpins->get_mesh();
  op_dat ls = mpins->get_ls();
  const DG_FP ls_alpha = mpins->get_ls_alpha();

  DGTempDat vol_0 = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat vol_1 = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat mass_0 = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat mass_1 = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat mass_total = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(measure_mass_of_phases, "measure_mass_of_phases", mesh->cells,
              op_arg_gbl(&ls_alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(ls, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mpins->rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vol_0.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vol_1.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(mass_0.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(mass_1.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(mass_total.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->mass(vol_0.dat);
  mesh->mass(vol_1.dat);
  mesh->mass(mass_0.dat);
  mesh->mass(mass_1.dat);
  mesh->mass(mass_total.dat);

  MassHistory data;
  data.time = mpins->get_time();
  data.phase0_vol = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&data.phase0_vol, 1, DG_FP_STR, OP_INC),
              op_arg_dat(vol_0.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  dg_dat_pool->releaseTempDatCells(vol_0);
  data.phase1_vol = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&data.phase1_vol, 1, DG_FP_STR, OP_INC),
              op_arg_dat(vol_1.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  dg_dat_pool->releaseTempDatCells(vol_1);
  data.phase0_mass = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&data.phase0_mass, 1, DG_FP_STR, OP_INC),
              op_arg_dat(mass_0.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  dg_dat_pool->releaseTempDatCells(mass_0);
  data.phase1_mass = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&data.phase1_mass, 1, DG_FP_STR, OP_INC),
              op_arg_dat(mass_1.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  dg_dat_pool->releaseTempDatCells(mass_1);
  data.total_mass = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&data.total_mass, 1, DG_FP_STR, OP_INC),
              op_arg_dat(mass_total.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  dg_dat_pool->releaseTempDatCells(mass_total);

  history.push_back(data);
}

std::string MassOfPhases2D::get_filename() {
  return "mass_of_phases";
}

std::string MassOfPhases2D::get_csv_header() {
  return "time,vol_phase_0,vol_phase_1,mass_phase_0,mass_phase_1,total_mass";
}

std::string MassOfPhases2D::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].phase0_vol) + ",";
    result = result + double_to_text(history[io_count].phase1_vol) + ",";
    result = result + double_to_text(history[io_count].phase0_mass) + ",";
    result = result + double_to_text(history[io_count].phase1_mass) + ",";
    result = result + double_to_text(history[io_count].total_mass);
    io_count++;
    return result;
  } else {
    return "";
  }
}