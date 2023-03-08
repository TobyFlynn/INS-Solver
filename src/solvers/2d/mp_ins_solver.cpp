#include "solvers/2d/mp_ins_solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"
#include "linear_solvers/petsc_amg.h"
#include "linear_solvers/petsc_block_jacobi.h"

extern Timing *timer;
extern DGConstants *constants;

using namespace std;

MPINSSolver2D::MPINSSolver2D(DGMesh2D *m) {
  // Hardcoded for the periodic cylinder case
  // int pressure_dirichlet[3] = {1, -1, -1};
  // int pressure_neumann[3] = {0, 2, -1};
  // int viscosity_dirichlet[3] = {0, 2, -1};
  // int viscosity_neumann[3] = {1, -1, -1};

  mesh = m;
  ls = new LevelSetSolver2D(mesh);
  pressureMatrix = new FactorPoissonMatrix2D(mesh);
  viscosityMatrix = new FactorMMPoissonMatrix2D(mesh);
  pressureSolver = new PETScAMGSolver(mesh);
  viscositySolver = new PETScBlockJacobiSolver(mesh);
  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_nullspace(true);
  viscositySolver->set_matrix(viscosityMatrix);
  viscositySolver->set_nullspace(false);

  // pressurePoisson->setDirichletBCs(pressure_dirichlet);
  // pressurePoisson->setNeumannBCs(pressure_neumann);
  // viscosityPoisson->setDirichletBCs(viscosity_dirichlet);
  // viscosityPoisson->setNeumannBCs(viscosity_neumann);

  std::string name;
  DG_FP *dg_np_data = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  for(int i = 0; i < 2; i++) {
    name = "mp_ins_solver_vel0" + std::to_string(i);
    vel[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "mp_ins_solver_vel1" + std::to_string(i);
    vel[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "mp_ins_solver_n0" + std::to_string(i);
    n[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "mp_ins_solver_n1" + std::to_string(i);
    n[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "mp_ins_solver_force0" + std::to_string(i);
    force[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "mp_ins_solver_force1" + std::to_string(i);
    force[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "mp_ins_solver_velT" + std::to_string(i);
    velT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "mp_ins_solver_velTT" + std::to_string(i);
    velTT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
  }
  for(int i = 0; i < 5; i++) {
    name = "mp_ins_solver_tmp_np" + std::to_string(i);
    tmp_np[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
  }
  pr  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, "mp_ins_solver_pr");
  rho = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, "mp_ins_solver_rho");
  mu  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, "mp_ins_solver_mu");
  free(dg_np_data);

  DG_FP *g_np_data = (DG_FP *)calloc(DG_G_NP * mesh->cells->size, sizeof(DG_FP));
  for(int i = 0; i < 5; i++) {
    string name    = "tmp_g_np" + to_string(i);
    tmp_g_np[i] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, g_np_data, name.c_str());
  }
  dPdN[0] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, g_np_data, "mp_ins_solver_dPdN0");
  dPdN[1] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, g_np_data, "mp_ins_solver_dPdN1");
  free(g_np_data);

  int *bc_1_data = (int *)calloc(mesh->bfaces->size, sizeof(int));
  bc_types     = op_decl_dat(mesh->bfaces, 1, "int", bc_1_data, "ins_solver_bc_types");
  pr_bc_types  = op_decl_dat(mesh->bfaces, 1, "int", bc_1_data, "ins_solver_pr_bc_types");
  vis_bc_types = op_decl_dat(mesh->bfaces, 1, "int", bc_1_data, "ins_solver_vis_bc_types");
  free(bc_1_data);

  DG_FP *tmp_np_np = (DG_FP *)calloc(DG_NP * DG_NP * mesh->cells->size, sizeof(DG_FP));
  proj_op_xx = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, tmp_np_np, "proj_op_xx");
  proj_op_yy = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, tmp_np_np, "proj_op_yy");
  proj_op_xy = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, tmp_np_np, "proj_op_xy");
  proj_op_yx = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, tmp_np_np, "proj_op_yx");
  free(tmp_np_np);

  DG_FP *tmp_DG_FP_1 = (DG_FP *)calloc(mesh->cells->size, sizeof(DG_FP));
  proj_pen = op_decl_dat(mesh->cells, 1, DG_FP_STR, tmp_DG_FP_1, "proj_pen");
  proj_h   = op_decl_dat(mesh->cells, 1, DG_FP_STR, tmp_DG_FP_1, "proj_h");
  free(tmp_DG_FP_1);

  f[0] = tmp_np[0];
  f[1] = tmp_np[1];
  f[2] = tmp_np[2];
  f[3] = tmp_np[3];
  divVelT = tmp_np[0];
  curlVel = tmp_np[1];
  gradCurlVel[0] = tmp_np[2];
  gradCurlVel[1] = tmp_np[3];
  pRHS = tmp_np[1];
  pr_mat_fact = tmp_np[2];
  dpdx = tmp_np[1];
  dpdy = tmp_np[2];
  proj_rhs_x = tmp_np[0];
  proj_rhs_y = tmp_np[3];
  visRHS[0] = tmp_np[0];
  visRHS[1] = tmp_np[1];
  vis_mat_mm_fact = tmp_np[2];
  ls_nx    = tmp_np[0];
  ls_ny    = tmp_np[1];
  ls_curv  = tmp_np[2];
  ls_delta = tmp_np[3];

  gVel[0] = tmp_g_np[0];
  gVel[1] = tmp_g_np[1];
  gAdvecFlux[0] = tmp_g_np[2];
  gAdvecFlux[1] = tmp_g_np[3];
  gN[0] = tmp_g_np[0];
  gN[1] = tmp_g_np[1];
  gGradCurl[0] = tmp_g_np[2];
  gGradCurl[1] = tmp_g_np[3];
  gRho = tmp_g_np[4];
  prBC = tmp_g_np[0];
  visBC[0] = tmp_g_np[0];
  visBC[1] = tmp_g_np[1];

  currentInd = 0;
  time = 0.0;

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;
}

MPINSSolver2D::~MPINSSolver2D() {
  delete pressureMatrix;
  delete viscosityMatrix;
  delete pressureSolver;
  delete viscositySolver;
  delete ls;
}

void MPINSSolver2D::init(const DG_FP re, const DG_FP refVel) {
  timer->startTimer("MPINSSolver2D - Init");
  reynolds = re;

  ls->init();

  // Set initial conditions
  op_par_loop(ins_2d_set_ic, "ins_2d_set_ic", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  dt = numeric_limits<DG_FP>::max();
  op_par_loop(calc_dt, "calc_dt", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_MIN));
  dt = dt / (DG_ORDER * DG_ORDER * refVel);
  op_printf("dt: %g\n", dt);

  if(mesh->bface2nodes) {
    op_par_loop(ins_bc_types, "ins_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,     -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(pr_bc_types,  -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

  ls->getRhoMu(rho, mu);
  ls->getNormalsCurvature(ls_nx, ls_ny, ls_curv);
  ls->getDiracDelta(ls_delta);

  op_par_loop(mp_ins_surf_ten_2d, "mp_ins_surf_ten_2d", mesh->cells,
              op_arg_dat(ls_nx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_ny, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_curv, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_delta, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(force[0][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(force[0][1], -1, OP_ID, DG_NP, "double", OP_WRITE));

  // Setup div-div pressure projection
  op_par_loop(project_2d_setup, "project_2d_setup", mesh->cells,
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(proj_op_xx, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(proj_op_yy, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(proj_op_yx, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(proj_op_xy, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE));

  op_par_loop(poisson_h, "poisson_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

  pressureMatrix->set_surface(ls->s, ls->alpha);
  viscosityMatrix->set_surface(ls->s, ls->alpha);

  timer->endTimer("MPINSSolver2D - Init");
}

void MPINSSolver2D::step() {
  timer->startTimer("MPINSSolver2D - Advection");
  advection();
  timer->endTimer("MPINSSolver2D - Advection");

  timer->startTimer("MPINSSolver2D - Pressure");
  pressure();
  timer->endTimer("MPINSSolver2D - Pressure");

  // timer->startTimer("Shock Capturing");
  // shock_capturing();
  // timer->endTimer("Shock Capturing");

  timer->startTimer("MPINSSolver2D - Viscosity");
  viscosity();
  timer->endTimer("MPINSSolver2D - Viscosity");

  // timer->startTimer("MPINSSolver2D - Surface");
  // surface();
  // timer->endTimer("MPINSSolver2D - Surface");

  currentInd = (currentInd + 1) % 2;
  time += dt;
  g0 = 1.5;
  a0 = 2.0;
  a1 = -0.5;
  b0 = 2.0;
  b1 = -1.0;
}

// Calculate Nonlinear Terms
void MPINSSolver2D::advection() {
  // Calculate flux values
  op_par_loop(ins_advec_flux_2d, "ins_advec_flux_2d", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0], f[1], n[currentInd][0]);
  mesh->div(f[2], f[3], n[currentInd][1]);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, vel[currentInd][0], 0.0, gVel[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, vel[currentInd][1], 0.0, gVel[1]);

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(gAdvecFlux[0], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(gAdvecFlux[1], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

  // Exchange values on edges between elements
  op_par_loop(ins_advec_faces_2d, "ins_advec_faces_2d", mesh->faces,
              op_arg_dat(mesh->order,     -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gVel[0],         -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gVel[1],         -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gAdvecFlux[0],   -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC),
              op_arg_dat(gAdvecFlux[1],   -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC));

  // Enforce BCs
  if(mesh->bface2cells) {
    op_par_loop(ins_advec_bc_2d, "ins_advec_bc_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->order,     0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->y,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gVel[0],         0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gVel[1],         0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gAdvecFlux[0],   0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC),
                op_arg_dat(gAdvecFlux[1],   0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC));
  }
  op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, gAdvecFlux[0], 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, gAdvecFlux[1], 1.0, n[currentInd][1]);

  // Calculate the intermediate velocity values
  op_par_loop(mp_ins_advec_intermediate_vel_2d, "mp_ins_advec_intermediate_vel_2d", mesh->cells,
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
              op_arg_dat(force[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(force[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}

bool MPINSSolver2D::pressure() {
  timer->startTimer("MPINSSolver2D - Pressure RHS");

  // mesh->div(velT[0], velT[1], divVelT);
  mesh->cub_div_with_central_flux(velT[0], velT[1], divVelT);
  mesh->curl(vel[currentInd][0], vel[currentInd][1], curlVel);
  mesh->grad(curlVel, gradCurlVel[0], gradCurlVel[1]);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, n[currentInd][0], 0.0, gN[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, n[currentInd][1], 0.0, gN[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, gradCurlVel[0], 0.0, gGradCurl[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, gradCurlVel[1], 0.0, gGradCurl[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, rho, 0.0, gRho);

  // Apply Neumann pressure boundary conditions
  if(mesh->bface2cells) {
    op_par_loop(mp_ins_pressure_bc_2d, "mp_ins_pressure_bc_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->y,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gN[0], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gN[1], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gGradCurl[0], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gGradCurl[1], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gRho, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dPdN[currentInd], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC));
  }
  // Apply Dirichlet BCs
  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(prBC, -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

  // Calculate RHS of pressure solve
  // This assumes that the boundaries will always be order DG_ORDER
  op_par_loop(ins_pressure_rhs_2d, "ins_pressure_rhs_2d", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->sJ, -1, OP_ID, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_RW),
              op_arg_dat(divVelT, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::MASS, divVelT, 0.0, pRHS);
  op2_gemv(mesh, true, 1.0, DGConstants::GAUSS_INTERP, dPdN[(currentInd + 1) % 2], 1.0, pRHS);
  timer->endTimer("MPINSSolver2D - Pressure RHS");

  // Call PETSc linear solver
  timer->startTimer("MPINSSolver2D - Pressure Linear Solve");

  op_par_loop(mp_ins_pressure_fact_2d, "mp_ins_pressure_fact_2d", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_mat_fact, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  bool converged;
  pressureMatrix->set_factor(pr_mat_fact);
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureMatrix->calc_mat();
  pressureSolver->set_bcs(prBC);
  converged = pressureSolver->solve(pRHS, pr);
  timer->endTimer("MPINSSolver2D - Pressure Linear Solve");

  timer->startTimer("MPINSSolver2D - Pressure Projection");

  project_velocity();

  timer->endTimer("MPINSSolver2D - Pressure Projection");
  return converged;
}

void MPINSSolver2D::project_velocity() {
  // Calculate gradient of pressure
  mesh->cub_grad_with_central_flux(pr, dpdx, dpdy);

  if(true) {
    // Calculate new velocity intermediate values
    op_par_loop(mp_ins_pressure_update_2d, "mp_ins_pressure_update_2d", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));
  } else {
    // Calculate new velocity intermediate values
    op_par_loop(mp_project_2d_0, "mp_project_2d_0", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(proj_rhs_x, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(proj_rhs_y, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

    mesh->mass(proj_rhs_x);
    mesh->mass(proj_rhs_y);

    DG_FP factor = dt * 1.0;
    // DG_FP factor = dt / Cr;
    // op_printf("Cr: %g\n", Cr);
    op_par_loop(project_2d_pen, "project_2d_pen", mesh->cells,
                op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

    // Do the 2 vector linear solve with project_mat as the matrix
    int num_cells = 0;
    int num_converge = 0;
    DG_FP num_iter = 0.0;
    op_par_loop(project_2d_cg, "project_2d_cg", mesh->cells,
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_op_xx, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_op_yy, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_op_yx, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_op_xy, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_rhs_x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_rhs_y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
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
  }
}

bool MPINSSolver2D::viscosity() {
  timer->startTimer("MPINSSolver2D - Viscosity RHS");
  DG_FP time_n1 = time + dt;

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(visBC[0], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visBC[1], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

  // Get BCs for viscosity solve
  if(mesh->bface2cells) {
    op_par_loop(ins_vis_bc_2d, "ins_vis_bc_2d", mesh->bfaces,
                op_arg_gbl(&time_n1, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->y,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(visBC[0], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC),
                op_arg_dat(visBC[1], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC));
  }
  // Set up RHS for viscosity solve
  op_par_loop(ins_vis_copy_2d, "ins_vis_copy_2d", mesh->cells,
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

   mesh->mass(visRHS[0]);
   mesh->mass(visRHS[1]);

  DG_FP factor = reynolds / dt;
  op_par_loop(mp_ins_vis_rhs_2d, "mp_ins_vis_rhs_2d", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(rho,       -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(visRHS[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  timer->endTimer("MPINSSolver2D - Viscosity RHS");

  // Call PETSc linear solver
  timer->startTimer("MPINSSolver2D - Viscosity Linear Solve");
  factor = g0 * reynolds / dt;
  op_par_loop(mp_ins_vis_mm_fact_2d, "mp_ins_vis_mm_fact_2d", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vis_mat_mm_fact, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  viscosityMatrix->set_factor(mu);
  viscosityMatrix->set_mm_factor(vis_mat_mm_fact);
  viscosityMatrix->set_bc_types(vis_bc_types);
  viscosityMatrix->calc_mat();
  viscositySolver->set_bcs(visBC[0]);
  bool convergedX = viscositySolver->solve(visRHS[0], vel[(currentInd + 1) % 2][0]);

  viscositySolver->set_bcs(visBC[1]);
  bool convergedY = viscositySolver->solve(visRHS[1], vel[(currentInd + 1) % 2][1]);
  timer->endTimer("MPINSSolver2D - Viscosity Linear Solve");

  // timer->startTimer("Filtering");
  // filter(mesh, data->Q[(currentInd + 1) % 2][0]);
  // filter(mesh, data->Q[(currentInd + 1) % 2][1]);
  // timer->endTimer("Filtering");

  return convergedX && convergedY;
}

void MPINSSolver2D::surface() {
  ls->setVelField(vel[(currentInd + 1) % 2][0], vel[(currentInd + 1) % 2][1]);
  ls->step(dt);
  ls->getRhoMu(rho, mu);
  ls->getNormalsCurvature(ls_nx, ls_ny, ls_curv);
  ls->getDiracDelta(ls_delta);

  op_par_loop(mp_ins_surf_ten_2d, "mp_ins_surf_ten_2d", mesh->cells,
              op_arg_dat(ls_nx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_ny, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_curv, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(ls_delta, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(force[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(force[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_WRITE));
}

// DG_FP MPINSSolver2D::getAvgPressureConvergance() {
//   return pressurePoisson->getAverageConvergeIter();
// }

// DG_FP MPINSSolver2D::getAvgViscosityConvergance() {
//   return viscosityPoisson->getAverageConvergeIter();
// }

DG_FP MPINSSolver2D::get_time() {
  return time;
}

DG_FP MPINSSolver2D::get_dt() {
  return dt;
}

void MPINSSolver2D::dump_data(const std::string &filename) {
  timer->startTimer("MPINSSolver2D - Dump Data");
  mesh->cub_div_with_central_flux(velT[0], velT[1], tmp_np[0]);
  mesh->cub_grad_with_central_flux(pr, tmp_np[1], tmp_np[2]);
  ls->getDiracDelta(tmp_np[3]);

  op_par_loop(fp_pen_term, "fp_pen_term", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&ls->alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 3 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(ls->s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_np[4], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(tmp_np[0], filename.c_str());
  op_fetch_data_hdf5_file(tmp_np[1], filename.c_str());
  op_fetch_data_hdf5_file(tmp_np[2], filename.c_str());
  op_fetch_data_hdf5_file(tmp_np[3], filename.c_str());
  op_fetch_data_hdf5_file(tmp_np[4], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][1], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][1], filename.c_str());
  op_fetch_data_hdf5_file(velT[0], filename.c_str());
  op_fetch_data_hdf5_file(velT[1], filename.c_str());
  op_fetch_data_hdf5_file(velTT[0], filename.c_str());
  op_fetch_data_hdf5_file(velTT[1], filename.c_str());
  op_fetch_data_hdf5_file(pr, filename.c_str());
  op_fetch_data_hdf5_file(rho, filename.c_str());
  op_fetch_data_hdf5_file(mu, filename.c_str());
  op_fetch_data_hdf5_file(ls->s, filename.c_str());
  op_fetch_data_hdf5_file(mesh->order, filename.c_str());
  timer->endTimer("MPINSSolver2D - Dump Data");
}
