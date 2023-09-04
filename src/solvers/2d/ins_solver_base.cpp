#include "solvers/2d/ins_solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "dg_utils.h"

#include "timing.h"
#include "config.h"
#include "dg_linear_solvers/petsc_amg.h"
#include "dg_linear_solvers/petsc_block_jacobi.h"
#include "dg_linear_solvers/petsc_pmultigrid.h"
#include "dg_linear_solvers/initial_guess_extrapolation.h"

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

INSSolverBase2D::INSSolverBase2D(DGMesh2D *m) {
  mesh = m;

  read_options();
  init_dats();

  std::string name;
  for(int i = 0; i < 2; i++) {
    name = "ins_solver_vel0" + std::to_string(i);
    vel[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_vel1" + std::to_string(i);
    vel[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_n0" + std::to_string(i);
    n[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_n1" + std::to_string(i);
    n[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
  }
  pr = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_pr");

  dPdN[0] = op_decl_dat(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_dPdN0");
  dPdN[1] = op_decl_dat(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_dPdN1");
}

INSSolverBase2D::INSSolverBase2D(DGMesh2D *m, const std::string &filename) {
  mesh = m;

  read_options();
  init_dats();

  vel[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel00");
  vel[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel10");
  vel[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel01");
  vel[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel11");
  n[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n00");
  n[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n10");
  n[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n01");
  n[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n11");
  pr = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_pr");
  dPdN[0] = op_decl_dat_hdf5(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, filename.c_str(), "ins_solver_dPdN0");
  dPdN[1] = op_decl_dat_hdf5(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, filename.c_str(), "ins_solver_dPdN1");
}

void INSSolverBase2D::read_options() {
  int tmp_div = 1;
  config->getInt("solver-options", "div_div", tmp_div);
  div_div_proj = tmp_div != 0;
  config->getInt("solver-options", "sub_cycle", sub_cycles);
  config->getInt("solver-options", "num_iter_before_sub_cycle", it_pre_sub_cycle);
  it_pre_sub_cycle = it_pre_sub_cycle > 1 ? it_pre_sub_cycle : 1;
  int tmp_eig = 1;
  config->getInt("solver-options", "extrapolate_initial_guess", tmp_eig);
  extrapolate_initial_guess = tmp_eig == 1;
  int tmp_shock = 0;
  config->getInt("solver-options", "shock_capturing", tmp_shock);
  shock_cap = tmp_shock == 1;

  filter_max_alpha = 18.0;
  config->getDouble("filter", "max_alpha", filter_max_alpha);
  filter_s0 = -2.0;
  config->getDouble("filter", "s0", filter_s0);
  filter_k = 1.0;
  config->getDouble("filter", "k", filter_k);
  filter_c = 1.0;
  config->getDouble("filter", "c", filter_c);

  int tmp_oia = 0;
  config->getInt("solver-options", "over_int_advec", tmp_oia);
  over_int_advec = tmp_oia == 1;
}

void INSSolverBase2D::init_dats() {
  std::string name;
  for(int i = 0; i < 2; i++) {
    name = "ins_solver_velT" + std::to_string(i);
    velT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_velTT" + std::to_string(i);
    velTT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
  }

  bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_bc_types");

  if(div_div_proj) {
    proj_h = op_decl_dat(mesh->cells, 1, DG_FP_STR, (DG_FP *)NULL, "ins_solver_proj_h");
  }
}

INSSolverBase2D::~INSSolverBase2D() {

}

void INSSolverBase2D::init(const DG_FP re, const DG_FP refVel) {
  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  op_printf("h: %g\n", h);

  if(div_div_proj) {
    DGTempDat tmp_npf = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
    op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
                op_arg_dat(tmp_npf.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

    op_par_loop(ins_proj_setup_0, "ins_proj_setup_0", mesh->faces,
                op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
                op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));

    op_par_loop(ins_proj_setup_1, "ins_proj_setup_1", mesh->cells,
                op_arg_dat(tmp_npf.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));
    dg_dat_pool->releaseTempDatCells(tmp_npf);
  }
}


void INSSolverBase2D::advec_current_non_linear() {
  DGTempDat f[2][2];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  // Calculate flux values
  op_par_loop(ins_advec_flux_2d, "ins_advec_flux_2d", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0][0].dat, f[0][1].dat, n[currentInd][0]);
  mesh->div(f[1][0].dat, f[1][1].dat, n[currentInd][1]);

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);

  DGTempDat tmp_advec_flux0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_advec_flux1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_2, "zero_npf_2", mesh->cells,
              op_arg_dat(tmp_advec_flux0.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(ins_advec_faces_2d, "ins_advec_faces_2d", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale,  -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0],  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1],  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_advec_flux0.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  // Enforce BCs
  if(mesh->bface2cells) {
    op_par_loop(ins_advec_bc_2d, "ins_advec_bc_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_advec_flux0.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_advec_flux1.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux0.dat, 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux1.dat, 1.0, n[currentInd][1]);

  dg_dat_pool->releaseTempDatCells(tmp_advec_flux0);
  dg_dat_pool->releaseTempDatCells(tmp_advec_flux1);
}

void INSSolverBase2D::advec_current_non_linear_over_int() {
  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int");
  DGTempDat f[2][2];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);

  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int - Interp");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, vel[currentInd][0], 0.0, f[0][0].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, vel[currentInd][1], 0.0, f[0][1].dat);
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int - Interp");

  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int - 0");
  op_par_loop(ins_2d_advec_oi_0, "ins_2d_advec_oi_0", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int - 0");

  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int - div");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, f[0][0].dat, 0.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, f[0][1].dat, 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, f[1][0].dat, 0.0, n[currentInd][1]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, f[1][1].dat, 1.0, n[currentInd][1]);
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int - div");

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);

  DGTempDat uM = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat vM = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat uP = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat vP = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(ins_2d_advec_oi_1, "ins_2d_advec_oi_1", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(vel[currentInd][0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(uM.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(vM.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(uP.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(vP.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  
  if(mesh->bface2cells) {
    op_par_loop(ins_2d_advec_oi_2, "ins_2d_advec_oi_2", mesh->bfaces,
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(vel[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(uM.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(vM.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(uP.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(vP.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));
  }

  DGTempDat uM_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat vM_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat uP_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat vP_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);

  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int - Interp Surf");
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, uM.dat, 0.0, uM_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, vM.dat, 0.0, vM_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, uP.dat, 0.0, uP_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, vP.dat, 0.0, vP_cub.dat);

  dg_dat_pool->releaseTempDatCells(uM);
  dg_dat_pool->releaseTempDatCells(vM);
  dg_dat_pool->releaseTempDatCells(uP);
  dg_dat_pool->releaseTempDatCells(vP);
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int - Interp Surf");

  op_par_loop(ins_2d_advec_oi_3, "ins_2d_advec_oi_3", mesh->cells,
              op_arg_dat(mesh->nx_c_new, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny_c_new, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ_c_new, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(uM_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vM_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(uP_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(vP_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, uP_cub.dat, -1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, vP_cub.dat, -1.0, n[currentInd][1]);

  dg_dat_pool->releaseTempDatCells(uM_cub);
  dg_dat_pool->releaseTempDatCells(vM_cub);
  dg_dat_pool->releaseTempDatCells(uP_cub);
  dg_dat_pool->releaseTempDatCells(vP_cub);
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int");
}

void INSSolverBase2D::advec_standard() {
  if(over_int_advec) {
    advec_current_non_linear_over_int();
  } else {
    advec_current_non_linear();
  }

  // Calculate the intermediate velocity values
  op_par_loop(ins_advec_intermediate_vel_2d, "ins_advec_intermediate_vel_2d", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}

void INSSolverBase2D::advec_standard(op_dat fx, op_dat fy, op_dat fx_old, op_dat fy_old) {
  if(over_int_advec) {
    advec_current_non_linear_over_int();
  } else {
    advec_current_non_linear();
  }

  // Calculate the intermediate velocity values
  op_par_loop(ins_advec_intermediate_vel_force_2d, "ins_advec_intermediate_vel_force_2d", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fx_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}

void INSSolverBase2D::advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat u_out,
                                      op_dat v_out, const double t) {
  double t0 = time - dt;
  double t1 = time;
  double tI = t;
  DGTempDat advec_sc[2];
  advec_sc[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_advec_sc_interp_2d, "ins_advec_sc_interp_2d", mesh->cells,
              op_arg_gbl(&t0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&t1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&tI, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  DGTempDat f[2][2];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_advec_sc_rhs_0_2d, "ins_advec_sc_rhs_0_2d", mesh->cells,
              op_arg_dat(u_in, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v_in, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0][0].dat, f[0][1].dat, u_out);
  mesh->div(f[1][0].dat, f[1][1].dat, v_out);

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);

  DGTempDat tmp_advec_flux0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_advec_flux1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_2, "zero_npf_2", mesh->cells,
              op_arg_dat(tmp_advec_flux0.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(ins_advec_sc_rhs_1_2d, "ins_advec_sc_rhs_1_2d", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale,  -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(u_in,  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v_in,  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[1].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_advec_flux0.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  if(mesh->bface2cells) {
    op_par_loop(ins_advec_sc_rhs_2_2d, "ins_advec_sc_rhs_2_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[0].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[1].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_advec_flux0.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_advec_flux1.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  dg_dat_pool->releaseTempDatCells(advec_sc[0]);
  dg_dat_pool->releaseTempDatCells(advec_sc[1]);

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux0.dat, 1.0, u_out);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux1.dat, 1.0, v_out);

  dg_dat_pool->releaseTempDatCells(tmp_advec_flux0);
  dg_dat_pool->releaseTempDatCells(tmp_advec_flux1);
}

void INSSolverBase2D::advec_sub_cycle_rk_step(const DG_FP time_sc, op_dat u, op_dat v) {
  // Request temporary dats
  DGTempDat advec_sc_rk[3][2];
  advec_sc_rk[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[2][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[2][1] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_advec_sc_copy_2d, "ins_advec_sc_copy_2d", mesh->cells,
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  for(int rk_step = 0; rk_step < 3; rk_step++) {
    double rk_time = time_sc;
    if(rk_step == 1) rk_time += sub_cycle_dt;
    if(rk_step == 2) rk_time += 0.5 * sub_cycle_dt;
    const int rk_ind = rk_step == 2 ? 2 : rk_step + 1;
    advec_sub_cycle_rhs(advec_sc_rk[0][0].dat, advec_sc_rk[0][1].dat,
                        advec_sc_rk[rk_ind][0].dat, advec_sc_rk[rk_ind][1].dat,
                        rk_time);
    // Set up next step
    if(rk_step == 0) {
      op_par_loop(ins_advec_sc_rk_0_2d, "ins_advec_sc_rk_0_2d", mesh->cells,
                  op_arg_gbl(&sub_cycle_dt, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                  op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    } else if(rk_step == 1) {
      op_par_loop(ins_advec_sc_rk_1_2d, "ins_advec_sc_rk_1_2d", mesh->cells,
                  op_arg_gbl(&sub_cycle_dt, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(advec_sc_rk[2][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[2][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                  op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    }
  }
  // Update velT
  op_par_loop(ins_advec_sc_rk_2_2d, "ins_advec_sc_rk_2_2d", mesh->cells,
              op_arg_gbl(&sub_cycle_dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[2][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[2][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  dg_dat_pool->releaseTempDatCells(advec_sc_rk[0][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[0][1]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[1][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[1][1]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[2][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[2][1]);
}

void INSSolverBase2D::advec_sub_cycle() {
  op_par_loop(ins_advec_sc_copy_2d, "ins_advec_sc_copy_2d", mesh->cells,
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // Advance 2 * number of subcycles
  for(int i = 0; i < 2 * sub_cycles; i++) {
    advec_sub_cycle_rk_step(time - dt + i * sub_cycle_dt, velT[0], velT[1]);
  }

  DGTempDat advec_sc_tmp[2];
  advec_sc_tmp[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_tmp[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  // Other velocity
  op_par_loop(ins_advec_sc_copy_2d, "ins_advec_sc_copy_2d", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // Advance number of subcycles
  for(int i = 0; i < sub_cycles; i++) {
    advec_sub_cycle_rk_step(time + i * sub_cycle_dt, advec_sc_tmp[0].dat, advec_sc_tmp[1].dat);
  }

  // Get final velT
  op_par_loop(ins_advec_sc_update_2d, "ins_advec_sc_update_2d", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  dg_dat_pool->releaseTempDatCells(advec_sc_tmp[0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_tmp[1]);

  // Needed for Pressure boundary conditions
  advec_current_non_linear();
}

void INSSolverBase2D::project_velocity(op_dat dpdx, op_dat dpdy) {
  if(!div_div_proj) {
    // Calculate new velocity intermediate values
    op_par_loop(ins_pressure_update_2d, "ins_pressure_update_2d", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  } else {
    DGTempDat projRHS[2];
    projRHS[0] = dg_dat_pool->requestTempDatCells(DG_NP);
    projRHS[1] = dg_dat_pool->requestTempDatCells(DG_NP);
    // Calculate new velocity intermediate values
    op_par_loop(ins_proj_rhs_2d, "ins_proj_rhs_2d", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    mesh->mass(projRHS[0].dat);
    mesh->mass(projRHS[1].dat);

    DGTempDat proj_pen = dg_dat_pool->requestTempDatCells(1);
    DG_FP factor = dt * 1.0;
    // DG_FP factor = dt / Cr;
    // op_printf("Cr: %g\n", Cr);
    op_par_loop(project_2d_pen, "project_2d_pen", mesh->cells,
                op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen.dat, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

    // Do the 2 vector linear solve with project_mat as the matrix
    int num_cells = 0;
    int num_converge = 0;
    DG_FP num_iter = 0.0;
    op_par_loop(ins_proj_cg_2d, "ins_proj_cg_2d", mesh->cells,
                op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen.dat, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(projRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(projRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_gbl(&num_cells, 1, "int", OP_INC),
                op_arg_gbl(&num_converge, 1, "int", OP_INC),
                op_arg_gbl(&num_iter, 1, DG_FP_STR, OP_INC));
    // op_printf("%d out of %d cells converged on projection step\n", num_converge, num_cells);
    // op_printf("Average iterations to converge on projection step %g\n", num_iter / (DG_FP)num_cells);
    if(num_cells != num_converge) {
      op_printf("%d out of %d cells converged on projection step\n", num_converge, num_cells);
      exit(-1);
    }

    dg_dat_pool->releaseTempDatCells(proj_pen);
    dg_dat_pool->releaseTempDatCells(projRHS[0]);
    dg_dat_pool->releaseTempDatCells(projRHS[1]);
  }
}

DG_FP INSSolverBase2D::max_vel() {
  DG_FP max_vel_tmp = 0.0;

  op_par_loop(ins_max_vel_2d, "ins_max_vel_2d", mesh->cells,
              op_arg_gbl(&max_vel_tmp, 1, DG_FP_STR, OP_MAX)
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  return std::max(1.0, sqrt(max_vel_tmp));
}

void INSSolverBase2D::add_to_pr_history() {
  const DG_FP t_n1 = time + dt;
  DGTempDat pr_copy = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(pr, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_copy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  pr_history.push_back({t_n1, pr_copy});

  while(pr_history.size() > EXTRAPOLATE_HISTORY_SIZE) {
    dg_dat_pool->releaseTempDatCells(pr_history[0].second);
    pr_history.erase(pr_history.begin());
  }
}

DG_FP INSSolverBase2D::get_time() {
  return time;
}

DG_FP INSSolverBase2D::get_dt() {
  return dt;
}

void INSSolverBase2D::shock_capture_filter_dat(op_dat in) {
  DGTempDat u_hat = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat u_modal = dg_dat_pool->requestTempDatCells(DG_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, in, 0.0, u_modal.dat);

  const double *r_ptr = constants->get_mat_ptr(DGConstants::R) + (DG_ORDER - 1) * DG_NP;
  const double *s_ptr = constants->get_mat_ptr(DGConstants::S) + (DG_ORDER - 1) * DG_NP;

  std::vector<DG_FP> r_vec, s_vec;
  for(int i = 0; i < DG_NP; i++) {
    r_vec.push_back(r_ptr[i]);
    s_vec.push_back(s_ptr[i]);
  }

  std::vector<DG_FP> simplex_vals = DGUtils::val_at_pt_N_1_2d_get_simplexes(r_vec, s_vec, DG_ORDER);

  op_par_loop(discont_sensor_0, "discont_sensor_0", mesh->cells,
              op_arg_gbl(simplex_vals.data(), DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(in, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_modal.dat, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_hat.dat, -1, OP_ID, DG_NP, "double", OP_WRITE));

  double max_alpha = filter_max_alpha;
  // double e0 = h;
  // double s0 = 1.0 / (double)(DG_ORDER * DG_ORDER * DG_ORDER * DG_ORDER);
  double s0 = filter_s0;
  // double k  = 5.0;
  double k = filter_k;
  double c = filter_c;
  op_par_loop(discont_sensor_cutoff_filter, "discont_sensor_cutoff_filter", mesh->cells,
              op_arg_gbl(&max_alpha, 1, "double", OP_READ),
              op_arg_gbl(&s0, 1, "double", OP_READ),
              op_arg_gbl(&k,  1, "double", OP_READ),
              op_arg_gbl(&c,  1, "double", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, "double", OP_READ),
              op_arg_dat(in, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_hat.dat, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_modal.dat, -1, OP_ID, DG_NP, "double", OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::V, u_modal.dat, 0.0, in);

  dg_dat_pool->releaseTempDatCells(u_hat);
  dg_dat_pool->releaseTempDatCells(u_modal);
}
