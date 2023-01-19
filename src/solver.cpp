#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

#include "solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_constants.h"
#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_operators.h"
#include "timing.h"

extern Timing *timer;
extern double dt, reynolds, refVel;

using namespace std;

void Solver::reverse_vel() {
  int tmp = 1;
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&tmp, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[0][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[0][1], -1, OP_ID, DG_NP, "double", OP_WRITE));
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&tmp, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[1][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[1][1], -1, OP_ID, DG_NP, "double", OP_WRITE));
}

Solver::Solver(std::string filename, int prob) {
  problem = prob;

  // Hardcoded for the periodic cylinder case
  int pressure_dirichlet[3] = {1, -1, -1};
  int pressure_neumann[3] = {0, 2, -1};
  int viscosity_dirichlet[3] = {0, 2, -1};
  int viscosity_neumann[3] = {1, -1, -1};

  mesh = new DGMesh(filename);

  data = new INSData(mesh);
  ls = new LS(mesh, data);
  pressurePoisson = new PetscPressureSolve(mesh, data, ls);
  viscosityPoisson = new PetscViscositySolve(mesh, data, ls);
  pMultigrid = new PMultigrid(mesh);

  pressurePoisson->setDirichletBCs(pressure_dirichlet);
  pressurePoisson->setNeumannBCs(pressure_neumann);
  viscosityPoisson->setDirichletBCs(viscosity_dirichlet);
  viscosityPoisson->setNeumannBCs(viscosity_neumann);
  pMultigrid->setDirichletBCs(pressure_dirichlet);
  pMultigrid->setNeumannBCs(pressure_neumann);

  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->edge2cells, NULL);

  timer->startTimer("Init functions");
  mesh->init();
  data->init();
  ls->init();
  pressurePoisson->init();
  viscosityPoisson->init();
  pMultigrid->init();
  timer->endTimer("Init functions");

  // Set initial conditions
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[0][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[0][1], -1, OP_ID, DG_NP, "double", OP_WRITE));
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[1][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[1][1], -1, OP_ID, DG_NP, "double", OP_WRITE));

  dt = numeric_limits<double>::max();
  op_par_loop(calc_dt, "calc_dt", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_MIN));
  dt = dt / (DG_ORDER * DG_ORDER * refVel);
  op_printf("dt: %g\n", dt);
}

Solver::~Solver() {
  delete pMultigrid;
  delete viscosityPoisson;
  delete pressurePoisson;
  delete ls;
  delete data;
  delete mesh;
}

void Solver::set_linear_solver(int ls) {
  linear_solver = ls;
}

// Calculate Nonlinear Terms
void Solver::advection(int currentInd, double a0, double a1, double b0,
                       double b1, double g0, double t) {
  // Calculate flux values
  op_par_loop(advection_flux, "advection_flux", mesh->cells,
              op_arg_dat(data->Q[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->F[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  div(mesh, data->F[0], data->F[1], data->N[currentInd][0]);
  div(mesh, data->F[2], data->F[3], data->N[currentInd][1]);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, data->Q[currentInd][0], 0.0, data->gQ[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, data->Q[currentInd][1], 0.0, data->gQ[1]);

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(data->flux[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(data->flux[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Exchange values on edges between elements
  op_par_loop(advection_faces, "advection_faces", mesh->edges,
              op_arg_dat(mesh->order,     -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(data->gQ[0],     -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[1],     -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->flux[0],   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->flux[1],   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  // Enforce BCs
  if(mesh->bedge2cells) {
    op_par_loop(advection_bc, "advection_bc", mesh->bedges,
                op_arg_dat(mesh->order,       0, mesh->bedge2cells, 1, "int", OP_READ),
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_gbl(&t, 1, "double", OP_READ),
                op_arg_gbl(&bc_time, 1, "double", OP_READ),
                op_arg_dat(mesh->gauss->x,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->y,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->gQ[0],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->gQ[1],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->flux[0],   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
                op_arg_dat(data->flux[1],   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));
  }
  op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[0], 1.0, data->N[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[1], 1.0, data->N[currentInd][1]);

  // Calculate the intermediate velocity values
  op_par_loop(advection_intermediate_vel, "advection_intermediate_vel", mesh->cells,
              op_arg_gbl(&a0, 1, "double", OP_READ),
              op_arg_gbl(&a1, 1, "double", OP_READ),
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->N[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->N[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->N[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->N[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QT[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->QT[1], -1, OP_ID, DG_NP, "double", OP_WRITE));
}

bool Solver::pressure(int currentInd, double a0, double a1, double b0,
                      double b1, double g0, double t) {
  timer->startTimer("Pressure Setup");

  div(mesh, data->QT[0], data->QT[1], data->divVelT);
  curl(mesh, data->Q[currentInd][0], data->Q[currentInd][1], data->curlVel);
  grad(mesh, data->curlVel, data->gradCurlVel[0], data->gradCurlVel[1]);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, data->N[currentInd][0], 0.0, data->gN[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, data->N[currentInd][1], 0.0, data->gN[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, data->gradCurlVel[0], 0.0, data->gGradCurl[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, data->gradCurlVel[1], 0.0, data->gGradCurl[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, data->rho, 0.0, data->gRho);

  // Apply pressure boundary conditions
  if(mesh->bedge2cells) {
    op_par_loop(pressure_bc, "pressure_bc", mesh->bedges,
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->order, 0, mesh->bedge2cells, 1, "int", OP_READ),
                op_arg_gbl(&t, 1, "double", OP_READ),
                op_arg_gbl(&bc_time, 1, "double", OP_READ),
                op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->gN[0], 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->gN[1], 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->gGradCurl[0], 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->gGradCurl[1], 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->gRho, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->dPdN[currentInd], 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));
  }

  // Calculate RHS of pressure solve
  // This assumes that the boundaries will always be order DG_ORDER
  op_par_loop(pressure_rhs, "pressure_rhs", mesh->cells,
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(data->divVelT, -1, OP_ID, DG_NP, "double", OP_RW));

  // op2_gemv(mesh, false, 1.0, DGConstants::LIFT, data->dPdN[(currentInd + 1) % 2], 1.0, data->divVelT);
  op2_gemv(mesh, false, 1.0, DGConstants::MASS, data->divVelT, 0.0, data->pRHS);
  op2_gemv(mesh, true, 1.0, DGConstants::GAUSS_INTERP, data->dPdN[(currentInd + 1) % 2], 1.0, data->pRHS);
  timer->endTimer("Pressure Setup");

  // Call PETSc linear solver
  timer->startTimer("Pressure Linear Solve");
  bool converged;
  if(linear_solver == 0) {
    pressurePoisson->setup();
    pressurePoisson->setBCValues(data->prBC);
    converged = pressurePoisson->solve(data->pRHS, data->p);
  } else {
    pMultigrid->setBCValues(data->prBC);
    pMultigrid->set_rho(data->rho);
    converged = pMultigrid->solve(data->pRHS, data->p);
  }
  timer->endTimer("Pressure Linear Solve");

  // Calculate gradient of pressure
  grad_with_central_flux(mesh, data->p, data->dpdx, data->dpdy);

  // Calculate new velocity intermediate values
  op_par_loop(pressure_update_vel, "pressure_update_vel", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->dpdx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->dpdy, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QT[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QT[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->QTT[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(data->prBC, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  return converged;
}

bool Solver::viscosity(int currentInd, double a0, double a1, double b0,
                       double b1, double g0, double t) {
  timer->startTimer("Viscosity Setup");
  double time = t + dt;

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(data->visBC[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(data->visBC[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Get BCs for viscosity solve
  if(mesh->bedge2cells) {
    op_par_loop(viscosity_bc, "viscosity_bc", mesh->bedges,
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_gbl(&time, 1, "double", OP_READ),
                op_arg_gbl(&bc_time, 1, "double", OP_READ),
                op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->visBC[0], 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
                op_arg_dat(data->visBC[1], 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));
  }
  // Set up RHS for viscosity solve
  /*
  op_par_loop(viscosity_mm, "viscosity_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->visRHS[0], -1, OP_ID, DG_NP, "double", OP_WRITE));
  op_par_loop(viscosity_mm, "viscosity_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->visRHS[1], -1, OP_ID, DG_NP, "double", OP_WRITE));
  */

  op_par_loop(viscosity_J, "viscosity_J", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->visTmp[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->visTmp[1], -1, OP_ID, DG_NP, "double", OP_WRITE));

   op2_gemv(mesh, false, 1.0, DGConstants::MASS, data->visTmp[0], 0.0, data->visRHS[0]);
   op2_gemv(mesh, false, 1.0, DGConstants::MASS, data->visTmp[1], 0.0, data->visRHS[1]);

  double factor = reynolds / dt;
  op_par_loop(viscosity_rhs, "viscosity_rhs", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(data->rho,       -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->visRHS[0], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(data->visRHS[1], -1, OP_ID, DG_NP, "double", OP_RW));

  timer->endTimer("Viscosity Setup");

  // Call PETSc linear solver
  timer->startTimer("Viscosity Linear Solve");
  factor = g0 * reynolds / dt;
  viscosityPoisson->setup(factor);
  viscosityPoisson->setBCValues(data->visBC[0]);
  bool convergedX = viscosityPoisson->solve(data->visRHS[0], data->Q[(currentInd + 1) % 2][0]);

  viscosityPoisson->setBCValues(data->visBC[1]);
  bool convergedY = viscosityPoisson->solve(data->visRHS[1], data->Q[(currentInd + 1) % 2][1]);
  timer->endTimer("Viscosity Linear Solve");

  // timer->startTimer("Filtering");
  // filter(mesh, data->Q[(currentInd + 1) % 2][0]);
  // filter(mesh, data->Q[(currentInd + 1) % 2][1]);
  // timer->endTimer("Filtering");

  return convergedX && convergedY;
}

void Solver::update_surface(int currentInd) {
  timer->startTimer("Surface");
  ls->setVelField(data->Q[(currentInd + 1) % 2][0], data->Q[(currentInd + 1) % 2][1]);
  ls->step(dt);
  timer->endTimer("Surface");
}

double Solver::getAvgPressureConvergance() {
  return pressurePoisson->getAverageConvergeIter();
}

double Solver::getAvgViscosityConvergance() {
  return viscosityPoisson->getAverageConvergeIter();
}

void Solver::set_bc_time(double t) {
  bc_time = t;
}

void Solver::switch_to_order(int o) {
  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);

  dats_to_update.push_back(data->N[0][0]);
  dats_to_update.push_back(data->N[0][1]);
  dats_to_update.push_back(data->N[1][0]);
  dats_to_update.push_back(data->N[1][1]);

  dats_to_update.push_back(data->p);

  dats_to_update.push_back(ls->s);

  dats_to_update.push_back(data->QT[0]);
  dats_to_update.push_back(data->QT[1]);
  dats_to_update.push_back(data->QTT[0]);
  dats_to_update.push_back(data->QTT[1]);

  mesh->update_order(o, dats_to_update);
  // ls->update_values();
}

void Solver::dump_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(data->Q[0][0], filename.c_str());
  op_fetch_data_hdf5_file(data->Q[0][1], filename.c_str());
  op_fetch_data_hdf5_file(data->Q[1][0], filename.c_str());
  op_fetch_data_hdf5_file(data->Q[1][1], filename.c_str());
  op_fetch_data_hdf5_file(data->QT[0], filename.c_str());
  op_fetch_data_hdf5_file(data->QT[1], filename.c_str());
  op_fetch_data_hdf5_file(data->QTT[0], filename.c_str());
  op_fetch_data_hdf5_file(data->QTT[1], filename.c_str());
  op_fetch_data_hdf5_file(data->p, filename.c_str());
  op_fetch_data_hdf5_file(data->rho, filename.c_str());
  op_fetch_data_hdf5_file(data->mu, filename.c_str());
  op_fetch_data_hdf5_file(ls->s, filename.c_str());
  op_fetch_data_hdf5_file(mesh->order, filename.c_str());
}