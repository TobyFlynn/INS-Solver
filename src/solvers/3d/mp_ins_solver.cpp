#include "solvers/3d/mp_ins_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_utils.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"
#include "utils.h"
#include "linear_solvers/petsc_amg.h"
#include "linear_solvers/petsc_block_jacobi.h"

#include <string>

extern Timing *timer;
extern DGConstants *constants;

MPINSSolver3D::MPINSSolver3D(DGMesh3D *m) {
  mesh = m;
  resuming = false;

  setup_common();

  std::string name;
  DG_FP * dg_np_data = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  for(int i = 0; i < 3; i++) {
    name = "ins_solver_vel0" + std::to_string(i);
    vel[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_vel1" + std::to_string(i);
    vel[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_n0" + std::to_string(i);
    n[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_n1" + std::to_string(i);
    n[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_force0" + std::to_string(i);
    force[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_force1" + std::to_string(i);
    force[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
  }
  pr  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, "ins_solver_pr");
  free(dg_np_data);

  DG_FP *dg_npf_data = (DG_FP *)calloc(4 * DG_NPF * mesh->cells->size, sizeof(DG_FP));
  dPdN[0]    = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_dPdN0");
  dPdN[1]    = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_dPdN1");
  free(dg_npf_data);

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  dt = 0.0;
  time = 0.0;

  currentInd = 0;

  lsSolver = new LevelSetSolver3D(mesh);
}

MPINSSolver3D::MPINSSolver3D(DGMesh3D *m, const std::string &filename, const int iter) {
  mesh = m;
  resuming = true;

  setup_common();

  vel[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel00");
  vel[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel10");
  vel[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel01");
  vel[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel11");
  vel[0][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel02");
  vel[1][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel12");
  n[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n00");
  n[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n10");
  n[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n01");
  n[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n11");
  n[0][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n02");
  n[1][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n12");
  force[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_force00");
  force[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_force10");
  force[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_force01");
  force[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_force11");
  force[0][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_force02");
  force[1][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_force12");
  pr = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_pr");
  dPdN[0] = op_decl_dat_hdf5(mesh->cells, 4 * DG_NPF, DG_FP_STR, filename.c_str(), "ins_solver_dPdN0");
  dPdN[1] = op_decl_dat_hdf5(mesh->cells, 4 * DG_NPF, DG_FP_STR, filename.c_str(), "ins_solver_dPdN1");

  currentInd = iter;

  if(iter > 0) {
    g0 = 1.5;
    a0 = 2.0;
    a1 = -0.5;
    b0 = 2.0;
    b1 = -1.0;
  } else {
    a0 = 1.0;
    a1 = 0.0;
    b0 = 1.0;
    b1 = 0.0;
    g0 = 1.0;
  }

  lsSolver = new LevelSetSolver3D(mesh, filename);
}

void MPINSSolver3D::setup_common() {
  std::string name;
  DG_FP * dg_np_data = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  for(int i = 0; i < 3; i++) {
    name = "ins_solver_velT" + std::to_string(i);
    velT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_velTT" + std::to_string(i);
    velTT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
  }
  for(int i = 0; i < 9; i++) {
    name = "ins_solver_tmp_np" + std::to_string(i);
    tmp_np[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
  }
  rho = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, "ins_solver_rho");
  mu  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, "ins_solver_mu");
  free(dg_np_data);

  DG_FP *dg_npf_data = (DG_FP *)calloc(4 * DG_NPF * mesh->cells->size, sizeof(DG_FP));
  tmp_npf[0] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_tmp_npf0");
  tmp_npf[1] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_tmp_npf1");
  tmp_npf[2] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_tmp_npf2");
  free(dg_npf_data);

  DG_FP *dg_npf_bc_data = (DG_FP *)calloc(DG_NPF * mesh->bfaces->size, sizeof(DG_FP));
  tmp_npf_bc = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, dg_npf_bc_data, "ins_solver_tmp_npf_bc");
  free(dg_npf_bc_data);

  int *bc_1_data = (int *)calloc(mesh->bfaces->size, sizeof(int));
  tmp_bc_1 = op_decl_dat(mesh->bfaces, 1, "int", bc_1_data, "ins_solver_tmp_bc_1");
  bc_types = op_decl_dat(mesh->bfaces, 1, "int", bc_1_data, "ins_solver_bc_types");
  free(bc_1_data);

  DG_FP *cell_1_data = (DG_FP *)calloc(mesh->cells->size, sizeof(DG_FP));
  art_vis = op_decl_dat(mesh->cells, 1, DG_FP_STR, cell_1_data, "ins_solver_art_vis");
  free(cell_1_data);

  f[0][0] = tmp_np[0]; f[0][1] = tmp_np[1]; f[0][2] = tmp_np[2];
  f[1][0] = tmp_np[3]; f[1][1] = tmp_np[4]; f[1][2] = tmp_np[5];
  f[2][0] = tmp_np[6]; f[2][1] = tmp_np[7]; f[2][2] = tmp_np[8];

  divVelT     = tmp_np[0];
  curlVel[0]  = tmp_np[1];
  curlVel[1]  = tmp_np[2];
  curlVel[2]  = tmp_np[3];
  curl2Vel[0] = tmp_np[4];
  curl2Vel[1] = tmp_np[5];
  curl2Vel[2] = tmp_np[6];
  pr_factor   = tmp_np[1];
  dpdx = tmp_np[0];
  dpdy = tmp_np[1];
  dpdz = tmp_np[2];
  shock_u = tmp_np[0];
  shock_u_hat = tmp_np[1];
  shock_u_modal = tmp_np[2];
  vis_factor    = tmp_np[0];
  vis_mm_factor = tmp_np[1];
  visRHS[0] = tmp_np[2];
  visRHS[1] = tmp_np[3];
  visRHS[2] = tmp_np[4];
  ls_normals[0] = tmp_np[0];
  ls_normals[1] = tmp_np[1];
  ls_normals[2] = tmp_np[2];
  ls_delta[0]   = tmp_np[3];
  ls_delta[1]   = tmp_np[4];
  ls_delta[2]   = tmp_np[5];
  ls_curv       = tmp_np[6];

  advec_flux[0] = tmp_npf[0];
  advec_flux[1] = tmp_npf[1];
  advec_flux[2] = tmp_npf[2];

  pr_bc_types  = tmp_bc_1;
  vis_bc_types = tmp_bc_1;
  pr_bc  = tmp_npf_bc;
  vis_bc = tmp_npf_bc;

  coarsePressureMatrix = new FactorPoissonCoarseMatrix3D(mesh);
  pressureMatrix = new FactorPoissonSemiMatrixFree3D(mesh);
  viscosityMatrix = new FactorMMPoissonSemiMatrixFree3D(mesh);
  // pressureSolver = new PETScAMGSolver(mesh);
  pressureSolver = new PETScPMultigrid(mesh);
  // pressureSolver = new PMultigridPoissonSolver(mesh);
  viscositySolver = new PETScBlockJacobiSolver(mesh);
  // viscositySolver = new PETScAMGSolver(mesh);

  pressureSolver->set_coarse_matrix(coarsePressureMatrix);
  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_nullspace(false);
  viscositySolver->set_matrix(viscosityMatrix);
  viscositySolver->set_nullspace(false);
}

MPINSSolver3D::~MPINSSolver3D() {
  delete coarsePressureMatrix;
  delete pressureMatrix;
  delete viscosityMatrix;
  delete pressureSolver;
  delete viscositySolver;
  delete lsSolver;
}

void MPINSSolver3D::init(const DG_FP re, const DG_FP refVel) {
  timer->startTimer("MPINSSolver3D - Init");
  reynolds = re;

  lsSolver->init();

  // Set initial conditions
  if(!resuming) {
    op_par_loop(ins_3d_set_ic, "ins_3d_set_ic", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1) * refVel);
  // dt *= 1e-2;
  op_printf("INS dt is %g\n", dt);
  time = dt * currentInd;
  currentInd = currentInd % 2;

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_bc_types, "ins_3d_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

  lsSolver->getRhoMu(rho, mu);
  lsSolver->getNormalsCurvature(ls_normals[0], ls_normals[1], ls_normals[2], ls_curv);
  lsSolver->getDiracDelta(ls_delta[0], ls_delta[1], ls_delta[2]);

  op_par_loop(mp_ins_3d_surf_ten, "mp_ins_3d_surf_ten", mesh->cells,
              op_arg_dat(ls_normals[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_normals[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_normals[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_curv, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_delta[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_delta[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_delta[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(force[0][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(force[0][1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(force[0][2], -1, OP_ID, DG_NP, "double", OP_WRITE));

  timer->endTimer("MPINSSolver3D - Init");
}

void MPINSSolver3D::step() {
  timer->startTimer("MPINSSolver3D - Advection");
  advection();
  timer->endTimer("MPINSSolver3D - Advection");

  timer->startTimer("MPINSSolver3D - Pressure");
  pressure();
  timer->endTimer("MPINSSolver3D - Pressure");

  // timer->startTimer("MPINSSolver3D - Shock Capturing");
  // shock_capturing();
  // timer->endTimer("MPINSSolver3D - Shock Capturing");

  timer->startTimer("MPINSSolver3D - Viscosity");
  viscosity();
  timer->endTimer("MPINSSolver3D - Viscosity");

  timer->startTimer("MPINSSolver3D - Surface");
  surface();
  timer->endTimer("MPINSSolver3D - Surface");

  currentInd = (currentInd + 1) % 2;
  time += dt;
  g0 = 1.5;
  a0 = 2.0;
  a1 = -0.5;
  b0 = 2.0;
  b1 = -1.0;
}

void MPINSSolver3D::advection() {
  op_par_loop(ins_3d_advec_0, "ins_3d_advec_0", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0][0], f[0][1], f[0][2], n[currentInd][0]);
  mesh->div(f[1][0], f[1][1], f[1][2], n[currentInd][1]);
  mesh->div(f[2][0], f[2][1], f[2][2], n[currentInd][2]);

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(advec_flux[0], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_flux[1], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_flux[2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

  // Flux across faces
  op_par_loop(ins_3d_advec_1, "ins_3d_advec_1", mesh->faces,
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_flux[0], -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(advec_flux[1], -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(advec_flux[2], -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));

  // Boundary flux
  if(mesh->bface2cells) {
    op_par_loop(ins_3d_advec_2, "ins_3d_advec_2", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_flux[0], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(advec_flux[1], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(advec_flux[2], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, advec_flux[0], 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, advec_flux[1], 1.0, n[currentInd][1]);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, advec_flux[2], 1.0, n[currentInd][2]);

  // Calculate the intermediate velocity values
  // TODO force
  op_par_loop(mp_ins_3d_advec_3, "mp_ins_3d_advec_3", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}

void MPINSSolver3D::pressure() {
  // mesh->div(velT[0], velT[1], velT[2], divVelT);
  mesh->div_with_central_flux(velT[0], velT[1], velT[2], divVelT);
  mesh->curl(vel[currentInd][0], vel[currentInd][1], vel[currentInd][2], curlVel[0], curlVel[1], curlVel[2]);
  // TODO potentially need to multiply curl by Mu here
  mesh->curl(curlVel[0], curlVel[1], curlVel[2], curl2Vel[0], curl2Vel[1], curl2Vel[2]);

  // TODO or maybe better to multiply by Mu after taking curl
  // (as how we smooth interface will have an effect before)
  if(mesh->bface2cells) {
    op_par_loop(mp_ins_3d_pr_0, "mp_ins_3d_pr_0", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dPdN[currentInd], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));
  }

  op_par_loop(mp_ins_3d_pr_1, "mp_ins_3d_pr_1", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(divVelT, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT);
  mesh->mass(divVelT);

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_pr_2, "ins_3d_pr_2", mesh->bfaces,
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(pr_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(pr_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  timer->startTimer("MPINSSolver3D - Pressure Linear Solve");
  pressureMatrix->set_factor(pr_factor);
  coarsePressureMatrix->set_factor(pr_factor);
  pressureMatrix->set_bc_types(pr_bc_types);
  coarsePressureMatrix->set_bc_types(pr_bc_types);
  pressureSolver->set_coarse_matrix(coarsePressureMatrix);
  pressureSolver->set_bcs(pr_bc);
  bool converged = pressureSolver->solve(divVelT, pr);
  if(!converged)
    throw std::runtime_error("\nPressure solve failed to converge\n");
  timer->endTimer("MPINSSolver3D - Pressure Linear Solve");

  mesh->grad_with_central_flux(pr, dpdx, dpdy, dpdz);

  op_par_loop(ins_3d_pr_3, "ins_3d_pr_3", mesh->cells,
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdz, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));
}

void MPINSSolver3D::viscosity() {
  DG_FP factor  = reynolds / dt;
  DG_FP factor2 = g0 * reynolds / dt;
  op_par_loop(mp_ins_3d_vis_0, "mp_ins_3d_vis_0", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&factor2, 1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mu,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(art_vis,  -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vis_factor,    -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vis_mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->mass(visRHS[0]);
  mesh->mass(visRHS[1]);
  mesh->mass(visRHS[2]);

  DG_FP vis_time = time + dt;
  if(mesh->bface2cells) {
    op_par_loop(ins_3d_vis_x, "ins_3d_vis_x", mesh->bfaces,
                op_arg_gbl(&vis_time, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  timer->startTimer("Vis Linear Solve");
  viscosityMatrix->set_factor(vis_factor);
  viscosityMatrix->set_mm_factor(vis_mm_factor);
  viscosityMatrix->set_bc_types(vis_bc_types);
  viscosityMatrix->calc_mat_partial();
  viscositySolver->set_bcs(vis_bc);
  bool convergedX = viscositySolver->solve(visRHS[0], vel[(currentInd + 1) % 2][0]);
  if(!convergedX)
    throw std::runtime_error("\nViscosity X solve failed to converge\n");

    if(mesh->bface2cells) {
      op_par_loop(ins_3d_vis_y, "ins_3d_vis_y", mesh->bfaces,
                  op_arg_gbl(&vis_time, 1, DG_FP_STR, OP_READ),
                  op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                  op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                  op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(velTT[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(velTT[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(velTT[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                  op_arg_dat(vis_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
    }

  viscositySolver->set_bcs(vis_bc);
  bool convergedY = viscositySolver->solve(visRHS[1], vel[(currentInd + 1) % 2][1]);
  if(!convergedY)
    throw std::runtime_error("\nViscosity Y solve failed to converge\n");

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_vis_z, "ins_3d_vis_z", mesh->bfaces,
                op_arg_gbl(&vis_time, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  viscositySolver->set_bcs(vis_bc);
  bool convergedZ = viscositySolver->solve(visRHS[2], vel[(currentInd + 1) % 2][2]);
  if(!convergedZ)
    throw std::runtime_error("\nViscosity Z solve failed to converge\n");
  timer->endTimer("Vis Linear Solve");
}

void MPINSSolver3D::surface() {
  lsSolver->setBCTypes(bc_types);
  // lsSolver->step(vel[(currentInd + 1) % 2][0], vel[(currentInd + 1) % 2][1], vel[(currentInd + 1) % 2][2], dt);
  lsSolver->getRhoMu(rho, mu);
  lsSolver->getNormalsCurvature(ls_normals[0], ls_normals[1], ls_normals[2], ls_curv);
  lsSolver->getDiracDelta(ls_delta[0], ls_delta[1], ls_delta[2]);

  op_par_loop(mp_ins_3d_surf_ten, "mp_ins_3d_surf_ten", mesh->cells,
              op_arg_dat(ls_normals[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_normals[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_normals[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_curv, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_delta[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_delta[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_delta[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(force[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(force[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(force[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, "double", OP_WRITE));
}

DG_FP MPINSSolver3D::get_time() {
  return time;
}

DG_FP MPINSSolver3D::get_dt() {
  return dt;
}

void MPINSSolver3D::dump_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(mesh->z, filename.c_str());
  op_fetch_data_hdf5_file(vel[0][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][1], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][2], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][1], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][2], filename.c_str());
  op_fetch_data_hdf5_file(n[0][0], filename.c_str());
  op_fetch_data_hdf5_file(n[0][1], filename.c_str());
  op_fetch_data_hdf5_file(n[0][2], filename.c_str());
  op_fetch_data_hdf5_file(n[1][0], filename.c_str());
  op_fetch_data_hdf5_file(n[1][1], filename.c_str());
  op_fetch_data_hdf5_file(n[1][2], filename.c_str());
  op_fetch_data_hdf5_file(velT[0], filename.c_str());
  op_fetch_data_hdf5_file(velT[1], filename.c_str());
  op_fetch_data_hdf5_file(velT[2], filename.c_str());
  op_fetch_data_hdf5_file(velTT[0], filename.c_str());
  op_fetch_data_hdf5_file(velTT[1], filename.c_str());
  op_fetch_data_hdf5_file(velTT[2], filename.c_str());
  op_fetch_data_hdf5_file(pr, filename.c_str());
  op_fetch_data_hdf5_file(rho, filename.c_str());
  op_fetch_data_hdf5_file(mu, filename.c_str());
  op_fetch_data_hdf5_file(dPdN[0], filename.c_str());
  op_fetch_data_hdf5_file(dPdN[1], filename.c_str());
  op_fetch_data_hdf5_file(lsSolver->s, filename.c_str());
}

void MPINSSolver3D::shock_capturing() {
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, pr, 0.0, shock_u_modal);

  const DG_FP *u_modal_ptr = getOP2PtrHost(shock_u_modal, OP_READ);
  const DG_FP *in_ptr = getOP2PtrHost(pr, OP_READ);
  DG_FP *u_ptr = getOP2PtrHost(shock_u, OP_WRITE);
  DG_FP *u_hat_ptr = getOP2PtrHost(shock_u_hat, OP_WRITE);

  // Get u_hat
  const DG_FP *r_ptr = constants->get_mat_ptr(DGConstants::R) + (DG_ORDER - 1) * constants->Np_max;
  const DG_FP *s_ptr = constants->get_mat_ptr(DGConstants::S) + (DG_ORDER - 1) * constants->Np_max;
  const DG_FP *t_ptr = constants->get_mat_ptr(DGConstants::T) + (DG_ORDER - 1) * constants->Np_max;
  #pragma omp parallel for
  for(int i = 0; i < mesh->cells->size; i++) {
    const DG_FP *modal = u_modal_ptr + i * DG_NP;
    for(int j = 0; j < DG_NP; j++) {
      u_ptr[i * DG_NP + j] = in_ptr[i * DG_NP + j];
      u_hat_ptr[i * DG_NP + j] = DGUtils::val_at_pt_N_1_3d(r_ptr[j], s_ptr[j], t_ptr[j], modal, DG_ORDER);
      u_hat_ptr[i * DG_NP + j] = u_ptr[i * DG_NP + j] - u_hat_ptr[i * DG_NP + j];
    }
  }

  releaseOP2PtrHost(shock_u_modal, OP_READ, u_modal_ptr);
  releaseOP2PtrHost(pr, OP_READ, in_ptr);
  releaseOP2PtrHost(shock_u, OP_WRITE, u_ptr);
  releaseOP2PtrHost(shock_u_hat, OP_WRITE, u_hat_ptr);


  // DG_FP e0 = h / (DG_FP)DG_ORDER;
  DG_FP e0 = h;
  DG_FP s0 = 1.0 / (DG_FP)(DG_ORDER * DG_ORDER * DG_ORDER * DG_ORDER);
  DG_FP k  = 5.0;
  // DG_FP k = 1.0;
  op_par_loop(discont_sensor, "discont_sensor", mesh->cells,
              op_arg_gbl(&e0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&s0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&k,  1, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(shock_u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(shock_u_hat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(art_vis, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));
}
