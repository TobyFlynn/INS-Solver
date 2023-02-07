#include "mp_ins_solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_op2_blas.h"

#include "timing.h"

extern Timing *timer;

using namespace std;

MPINSSolver2D::MPINSSolver2D(DGMesh2D *m) {
  // Hardcoded for the periodic cylinder case
  int pressure_dirichlet[3] = {1, -1, -1};
  int pressure_neumann[3] = {0, 2, -1};
  int viscosity_dirichlet[3] = {0, 2, -1};
  int viscosity_neumann[3] = {1, -1, -1};

  mesh = m;
  ls = new LevelSetSolver2D(mesh);
  pressurePoisson = new PetscPressureSolve(mesh);
  viscosityPoisson = new PetscViscositySolve(mesh);

  pressurePoisson->setDirichletBCs(pressure_dirichlet);
  pressurePoisson->setNeumannBCs(pressure_neumann);
  viscosityPoisson->setDirichletBCs(viscosity_dirichlet);
  viscosityPoisson->setNeumannBCs(viscosity_neumann);

  std::string name;
  double *dg_np_data = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  for(int i = 0; i < 2; i++) {
    name = "mp_ins_solver_vel0" + std::to_string(i);
    vel[0][i] = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, name.c_str());
    name = "mp_ins_solver_vel1" + std::to_string(i);
    vel[1][i] = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, name.c_str());
    name = "mp_ins_solver_n0" + std::to_string(i);
    n[0][i] = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, name.c_str());
    name = "mp_ins_solver_n1" + std::to_string(i);
    n[1][i] = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, name.c_str());
    name = "mp_ins_solver_velT" + std::to_string(i);
    velT[i] = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, name.c_str());
    name = "mp_ins_solver_velTT" + std::to_string(i);
    velTT[i] = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, name.c_str());
  }
  for(int i = 0; i < 4; i++) {
    name = "mp_ins_solver_tmp_np" + std::to_string(i);
    tmp_np[i] = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, name.c_str());
  }
  pr  = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, "mp_ins_solver_pr");
  rho = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, "mp_ins_solver_rho");
  mu  = op_decl_dat(mesh->cells, DG_NP, "double", dg_np_data, "mp_ins_solver_mu");
  free(dg_np_data);

  double *g_np_data = (double *)calloc(DG_G_NP * mesh->cells->size, sizeof(double));
  for(int i = 0; i < 5; i++) {
    string name    = "tmp_g_np" + to_string(i);
    tmp_g_np[i] = op_decl_dat(mesh->cells, DG_G_NP, "double", g_np_data, name.c_str());
  }
  dPdN[0] = op_decl_dat(mesh->cells, DG_G_NP, "double", g_np_data, "mp_ins_solver_dPdN0");
  dPdN[1] = op_decl_dat(mesh->cells, DG_G_NP, "double", g_np_data, "mp_ins_solver_dPdN1");
  free(g_np_data);

  f[0] = tmp_np[0];
  f[1] = tmp_np[1];
  f[2] = tmp_np[2];
  f[3] = tmp_np[3];
  divVelT = tmp_np[0];
  curlVel = tmp_np[1];
  gradCurlVel[0] = tmp_np[2];
  gradCurlVel[1] = tmp_np[3];
  pRHS = tmp_np[1];
  dpdx = tmp_np[1];
  dpdy = tmp_np[2];
  visRHS[0] = tmp_np[0];
  visRHS[1] = tmp_np[1];

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
  delete viscosityPoisson;
  delete pressurePoisson;
  delete ls;
}

void MPINSSolver2D::init(const double re, const double refVel) {
  timer->startTimer("MPINS - Init");
  reynolds = re;

  ls->init();
  pressurePoisson->init();
  viscosityPoisson->init();

  // Set initial conditions
  op_par_loop(ins_set_ic, "ins_set_ic", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(vel[0][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(vel[0][1], -1, OP_ID, DG_NP, "double", OP_WRITE));

  op_par_loop(ins_set_ic, "ins_set_ic", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(vel[1][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(vel[1][1], -1, OP_ID, DG_NP, "double", OP_WRITE));

  dt = numeric_limits<double>::max();
  op_par_loop(calc_dt, "calc_dt", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_MIN));
  dt = dt / (DG_ORDER * DG_ORDER * refVel);
  op_printf("dt: %g\n", dt);

  ls->getRhoMu(rho, mu);

  timer->endTimer("MPINS - Init");
}

void MPINSSolver2D::step() {
  timer->startTimer("Advection");
  advection();
  timer->endTimer("Advection");

  timer->startTimer("Pressure");
  pressure();
  timer->endTimer("Pressure");

  // timer->startTimer("Shock Capturing");
  // shock_capturing();
  // timer->endTimer("Shock Capturing");

  timer->startTimer("Viscosity");
  viscosity();
  timer->endTimer("Viscosity");

  timer->startTimer("Surface");
  surface();
  timer->endTimer("Surface");

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
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(f[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(f[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(f[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(f[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  mesh->div(f[0], f[1], n[currentInd][0]);
  mesh->div(f[2], f[3], n[currentInd][1]);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, vel[currentInd][0], 0.0, gVel[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, vel[currentInd][1], 0.0, gVel[1]);

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(gAdvecFlux[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(gAdvecFlux[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Exchange values on edges between elements
  op_par_loop(ins_advec_faces_2d, "ins_advec_faces_2d", mesh->faces,
              op_arg_dat(mesh->order,     -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gVel[0],         -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gVel[1],         -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gAdvecFlux[0],   -2, mesh->face2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(gAdvecFlux[1],   -2, mesh->face2cells, DG_G_NP, "double", OP_INC));

  // Enforce BCs
  if(mesh->bface2cells) {
    op_par_loop(ins_advec_bc_2d, "ins_advec_bc_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, "double", OP_READ),
                op_arg_dat(mesh->order,       0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x,  0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->y,  0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(gVel[0],         0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(gVel[1],         0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(gAdvecFlux[0],   0, mesh->bface2cells, DG_G_NP, "double", OP_INC),
                op_arg_dat(gAdvecFlux[1],   0, mesh->bface2cells, DG_G_NP, "double", OP_INC));
  }
  op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, gAdvecFlux[0], 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, gAdvecFlux[1], 1.0, n[currentInd][1]);

  // Calculate the intermediate velocity values
  op_par_loop(ins_advec_intermediate_vel_2d, "ins_advec_intermediate_vel_2d", mesh->cells,
              op_arg_gbl(&a0, 1, "double", OP_READ),
              op_arg_gbl(&a1, 1, "double", OP_READ),
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, "double", OP_WRITE));
}

bool MPINSSolver2D::pressure() {
  timer->startTimer("Pressure Setup");

  mesh->div(velT[0], velT[1], divVelT);
  mesh->curl(vel[currentInd][0], vel[currentInd][1], curlVel);
  mesh->grad(curlVel, gradCurlVel[0], gradCurlVel[1]);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, n[currentInd][0], 0.0, gN[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, n[currentInd][1], 0.0, gN[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, gradCurlVel[0], 0.0, gGradCurl[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, gradCurlVel[1], 0.0, gGradCurl[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, rho, 0.0, gRho);

  // Apply Neumann pressure boundary conditions
  if(mesh->bface2cells) {
    op_par_loop(ins_pressure_bc_2d, "ins_pressure_bc_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, "double", OP_READ),
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x,  0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->y,  0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(gN[0], 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(gN[1], 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(gGradCurl[0], 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(gGradCurl[1], 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(gRho, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(dPdN[currentInd], 0, mesh->bface2cells, DG_G_NP, "double", OP_INC));
  }
  // Apply Dirichlet BCs
  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(prBC, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Calculate RHS of pressure solve
  // This assumes that the boundaries will always be order DG_ORDER
  op_par_loop(ins_pressure_rhs_2d, "ins_pressure_rhs_2d", mesh->cells,
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(divVelT, -1, OP_ID, DG_NP, "double", OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::MASS, divVelT, 0.0, pRHS);
  op2_gemv(mesh, true, 1.0, DGConstants::GAUSS_INTERP, dPdN[(currentInd + 1) % 2], 1.0, pRHS);
  timer->endTimer("Pressure Setup");

  // Call PETSc linear solver
  timer->startTimer("Pressure Linear Solve");
  bool converged;
  pressurePoisson->setup(rho);
  pressurePoisson->setBCValues(prBC);
  converged = pressurePoisson->solve(pRHS, pr);
  timer->endTimer("Pressure Linear Solve");

  // Calculate gradient of pressure
  mesh->grad_with_central_flux(pr, dpdx, dpdy);

  // Calculate new velocity intermediate values
  op_par_loop(ins_pressure_update_2d, "ins_pressure_update_2d", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dpdx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dpdy, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  return converged;
}

bool MPINSSolver2D::viscosity() {
  timer->startTimer("Viscosity Setup");
  double time_n1 = time + dt;

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(visBC[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(visBC[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Get BCs for viscosity solve
  if(mesh->bface2cells) {
    op_par_loop(ins_vis_bc_2d, "ins_vis_bc_2d", mesh->bfaces,
                op_arg_gbl(&time_n1, 1, "double", OP_READ),
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->y, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(visBC[0], 0, mesh->bface2cells, DG_G_NP, "double", OP_INC),
                op_arg_dat(visBC[1], 0, mesh->bface2cells, DG_G_NP, "double", OP_INC));
  }
  // Set up RHS for viscosity solve
  op_par_loop(ins_vis_copy_2d, "ins_vis_copy_2d", mesh->cells,
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(visRHS[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(visRHS[1], -1, OP_ID, DG_NP, "double", OP_WRITE));

   mesh->mass(visRHS[0]);
   mesh->mass(visRHS[1]);

  double factor = reynolds / dt;
  op_par_loop(ins_vis_rhs_2d, "ins_vis_rhs_2d", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(rho,       -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(visRHS[0], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(visRHS[1], -1, OP_ID, DG_NP, "double", OP_RW));

  timer->endTimer("Viscosity Setup");

  // Call PETSc linear solver
  timer->startTimer("Viscosity Linear Solve");
  factor = g0 * reynolds / dt;
  viscosityPoisson->setup(factor, rho, mu);
  viscosityPoisson->setBCValues(visBC[0]);
  bool convergedX = viscosityPoisson->solve(visRHS[0], vel[(currentInd + 1) % 2][0]);

  viscosityPoisson->setBCValues(visBC[1]);
  bool convergedY = viscosityPoisson->solve(visRHS[1], vel[(currentInd + 1) % 2][1]);
  timer->endTimer("Viscosity Linear Solve");

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
}

double MPINSSolver2D::getAvgPressureConvergance() {
  return pressurePoisson->getAverageConvergeIter();
}

double MPINSSolver2D::getAvgViscosityConvergance() {
  return viscosityPoisson->getAverageConvergeIter();
}

void MPINSSolver2D::dump_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
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
}
