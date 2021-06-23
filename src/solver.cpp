#include "solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_constants.h"
#include "dg_blas_calls.h"
#include "dg_operators.h"
#include "load_mesh.h"

extern Timing *timer;
extern DGConstants *constants;

using namespace std;

Solver::Solver(std::string filename, int pmethod, int prob, bool multi) {
  problem = prob;
  multiphase = multi;

  // Ownership of the pointers is passed to DGMesh
  // so don't have to worry about freeing them
  double *coords_data;
  int *cells_data, *edge2node_data, *edge2cell_data, *bedge2node_data;
  int *bedge2cell_data, *bedge_type_data, *edgeNum_data, *bedgeNum_data;
  int numNodes_g, numCells_g, numEdges_g, numBoundaryEdges_g, numNodes;
  int numCells, numEdges, numBoundaryEdges;

  int pressure_dirichlet[3];
  int pressure_neumann[3];
  int viscosity_dirichlet[3];
  int viscosity_neumann[3];

  load_mesh(filename, &coords_data, &cells_data, &edge2node_data,
            &edge2cell_data, &bedge2node_data, &bedge2cell_data,
            &bedge_type_data, &edgeNum_data, &bedgeNum_data, &numNodes_g,
            &numCells_g, &numEdges_g, &numBoundaryEdges_g, &numNodes, &numCells,
            &numEdges, &numBoundaryEdges, pressure_dirichlet, pressure_neumann,
            viscosity_dirichlet, viscosity_neumann);

  mesh = new DGMesh(coords_data, cells_data, edge2node_data, edge2cell_data,
                    bedge2node_data, bedge2cell_data, bedge_type_data,
                    edgeNum_data, bedgeNum_data, numNodes_g, numCells_g,
                    numEdges_g, numBoundaryEdges_g, numNodes, numCells,
                    numEdges, numBoundaryEdges);
  data = new INSData(mesh);
  if(multiphase) {
    ls = new LS(mesh, data);
  }

  if(pmethod == 0) {
    Poisson_M *pressureM = new Poisson_M(mesh, data);
    pressureM->setDirichletBCs(pressure_dirichlet);
    pressureM->setNeumannBCs(pressure_neumann);
    pressurePoisson = pressureM;
    Poisson_M *viscosityM = new Poisson_M(mesh, data);
    viscosityM->setDirichletBCs(viscosity_dirichlet);
    viscosityM->setNeumannBCs(viscosity_neumann);
    viscosityPoisson = viscosityM;
  } else if(pmethod == 1) {
    Poisson_MF *pressureMF = new Poisson_MF(mesh, data);
    pressureMF->setDirichletBCs(pressure_dirichlet);
    pressureMF->setNeumannBCs(pressure_neumann);
    pressurePoisson = pressureMF;
    Poisson_MF *viscosityMF = new Poisson_MF(mesh, data);
    viscosityMF->setDirichletBCs(viscosity_dirichlet);
    viscosityMF->setNeumannBCs(viscosity_neumann);
    viscosityPoisson = viscosityMF;
  } else {
    Poisson_MF2 *pressureMF2 = new Poisson_MF2(mesh, data);
    pressureMF2->setDirichletBCs(pressure_dirichlet);
    pressureMF2->setNeumannBCs(pressure_neumann);
    pressurePoisson = pressureMF2;
    Poisson_MF2 *viscosityMF2 = new Poisson_MF2(mesh, data);
    viscosityMF2->setDirichletBCs(viscosity_dirichlet);
    viscosityMF2->setNeumannBCs(viscosity_neumann);
    viscosityPoisson = viscosityMF2;
  }

  op_partition("PARMETIS", "KWAY", mesh->cells, mesh->edge2cells, NULL);

  mesh->init();
  data->init();
  if(multiphase) {
    ls->init();
  }
  pressurePoisson->init();
  viscosityPoisson->init();

  // Set initial conditions
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x,   -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(mesh->y,   -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->nu,  -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[0][0],   -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->Q[0][1],   -1, OP_ID, 15, "double", OP_WRITE));

  dt = numeric_limits<double>::max();
  op_par_loop(calc_dt, "calc_dt", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_MIN));
  dt = dt / 25.0;
  op_printf("dt: %g\n", dt);
}

Solver::~Solver() {
  delete viscosityPoisson;
  delete pressurePoisson;
  if(multiphase) {
    delete ls;
  }
  delete data;
  delete mesh;
}

void Solver::advection(int currentInd, double a0, double a1, double b0,
                       double b1, double g0, double t) {
  // Calculate flux values
  op_par_loop(advection_flux, "advection_flux", mesh->cells,
              op_arg_dat(data->Q[currentInd][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->F[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[3], -1, OP_ID, 15, "double", OP_WRITE));

  div(mesh, data->F[0], data->F[1], data->N[currentInd][0]);
  div(mesh, data->F[2], data->F[3], data->N[currentInd][1]);

  // Exchange values on edges between elements
  op_par_loop(advection_faces, "advection_faces", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(data->Q[currentInd][0], -2, mesh->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -2, mesh->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->exQ[0], -2, mesh->edge2cells, 15, "double", OP_INC),
              op_arg_dat(data->exQ[1], -2, mesh->edge2cells, 15, "double", OP_INC));

  // Enforce BCs
  op_par_loop(advection_bc, "advection_bc", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->nu, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][0], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->exQ[0], 0, mesh->bedge2cells, 15, "double", OP_INC),
              op_arg_dat(data->exQ[1], 0, mesh->bedge2cells, 15, "double", OP_INC));

  // Calculate numberical flux across edges
  op_par_loop(advection_numerical_flux, "advection_numerical_flux", mesh->cells,
              op_arg_dat(mesh->fscale,  -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(mesh->nx,      -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(mesh->ny,      -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->exQ[0],  -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->exQ[1],  -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->flux[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->flux[1], -1, OP_ID, 15, "double", OP_WRITE));

  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(DGConstants::LIFT), 15, data->flux[0], 1.0, data->N[currentInd][0]);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(DGConstants::LIFT), 15, data->flux[1], 1.0, data->N[currentInd][1]);

  // Calculate the intermediate velocity values
  op_par_loop(advection_intermediate_vel, "advection_intermediate_vel", mesh->cells,
              op_arg_gbl(&a0, 1, "double", OP_READ),
              op_arg_gbl(&a1, 1, "double", OP_READ),
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&g0, 1, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[(currentInd + 1) % 2][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[(currentInd + 1) % 2][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->N[currentInd][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->N[currentInd][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->N[(currentInd + 1) % 2][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->N[(currentInd + 1) % 2][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->QT[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->QT[1], -1, OP_ID, 15, "double", OP_WRITE));
}

bool Solver::pressure(int currentInd, double a0, double a1, double b0,
                      double b1, double g0, double t) {
  timer->startPressureSetup();
  div(mesh, data->QT[0], data->QT[1], data->divVelT);
  curl(mesh, data->Q[currentInd][0], data->Q[currentInd][1], data->curlVel);
  grad(mesh, data->curlVel, data->gradCurlVel[0], data->gradCurlVel[1]);

  // Apply pressure boundary conditions
  op_par_loop(pressure_bc, "pressure_bc", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(mesh->nx, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(mesh->ny, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->nu, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->N[currentInd][0], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->N[currentInd][1], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->gradCurlVel[0], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->gradCurlVel[1], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], 0, mesh->bedge2cells, 15, "double", OP_INC));

  if(problem == 1) {
    op_par_loop(pressure_bc2, "pressure_bc2", mesh->bedges,
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_gbl(&t, 1, "double", OP_READ),
                op_arg_gbl(&problem, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, 21, "double", OP_READ),
                op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, 21, "double", OP_READ),
                op_arg_dat(data->gNu, 0, mesh->bedge2cells, 21, "double", OP_READ),
                op_arg_dat(data->prBC, 0, mesh->bedge2cells, 21, "double", OP_INC));
  }

  // Calculate RHS of pressure solve
  op_par_loop(pressure_rhs, "pressure_rhs", mesh->cells,
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&g0, 1, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->divVelT, -1, OP_ID, 15, "double", OP_RW));

  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(DGConstants::LIFT), 15, data->dPdN[(currentInd + 1) % 2], 1.0, data->divVelT);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(DGConstants::MASS), 15, data->divVelT, 0.0, data->pRHS);
  timer->endPressureSetup();

  // Call PETSc linear solver
  timer->startPressureLinearSolve();
  pressurePoisson->setBCValues(data->prBC);
  bool converged = pressurePoisson->solve(data->pRHS, data->p);
  timer->endPressureLinearSolve();

  // Calculate gradient of pressure
  grad(mesh, data->p, data->dpdx, data->dpdy);

  // Calculate new velocity intermediate values
  double factor = dt / g0;
  op_par_loop(pressure_update_vel, "pressure_update_vel", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->dpdx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->dpdy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->QT[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->QT[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->QTT[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->QTT[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->prBC, -1, OP_ID, 21, "double", OP_WRITE));

  return converged;
}

bool Solver::viscosity(int currentInd, double a0, double a1, double b0,
                       double b1, double g0, double t) {
  timer->startViscositySetup();
  double time = t + dt;
  // Get BCs for viscosity solve
  op_par_loop(viscosity_bc, "viscosity_bc", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&time, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(data->gNu, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(data->visBC[0], 0, mesh->bedge2cells, 21, "double", OP_INC),
              op_arg_dat(data->visBC[1], 0, mesh->bedge2cells, 21, "double", OP_INC));

  // Set up RHS for viscosity solve
  op2_gemv_batch(false, 15, 15, 1.0, mesh->cubature->mm, 15, data->QTT[0], 0.0, data->visRHS[0]);
  op2_gemv_batch(false, 15, 15, 1.0, mesh->cubature->mm, 15, data->QTT[1], 0.0, data->visRHS[1]);

  double factor;
  if(multiphase) {
    // factor = g0 / dt;
    factor = g0 / (nu0 * dt);
  } else {
    factor = g0 / (nu0 * dt);
  }

  op_par_loop(viscosity_rhs, "viscosity_rhs", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(data->visRHS[0], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->visRHS[1], -1, OP_ID, 15, "double", OP_RW));

  if(multiphase) {
    op_par_loop(viscosity_rhs_rho, "viscosity_rhs_rho", mesh->cells,
                op_arg_dat(data->rho, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(data->nu, -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(data->visRHS[0], -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(data->visRHS[1], -1, OP_ID, 15, "double", OP_RW));
  }

  timer->endViscositySetup();

  // Call PETSc linear solver
  timer->startViscosityLinearSolve();
  viscosityPoisson->setBCValues(data->visBC[0]);
  bool convergedX = viscosityPoisson->solve(data->visRHS[0], data->Q[(currentInd + 1) % 2][0], true, factor);

  viscosityPoisson->setBCValues(data->visBC[1]);
  bool convergedY = viscosityPoisson->solve(data->visRHS[1], data->Q[(currentInd + 1) % 2][1], true, factor);
  timer->endViscosityLinearSolve();

  // Reset BC dats ready for next iteration
  op_par_loop(viscosity_reset_bc, "viscosity_reset_bc", mesh->cells,
              op_arg_dat(data->visBC[0], -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(data->visBC[1], -1, OP_ID, 21, "double", OP_WRITE));

  return convergedX && convergedY;
}

void Solver::update_surface(int currentInd) {
  timer->startSurface();
  if(multiphase) {
    ls->setVelField(data->Q[(currentInd + 1) % 2][0], data->Q[(currentInd + 1) % 2][1]);
    ls->step(dt);
  }
  timer->endSurface();
}

// Function to calculate lift and drag coefficients of the cylinder
void Solver::lift_drag_coeff(double *lift, double *drag, int ind) {
  *lift = 0.0;
  *drag = 0.0;

  grad(mesh, data->Q[(ind + 1) % 2][0], data->dQdx[0], data->dQdy[0]);
  grad(mesh, data->Q[(ind + 1) % 2][1], data->dQdx[1], data->dQdy[1]);

  op_par_loop(lift_drag, "lift_drag", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->p, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dQdx[0], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dQdy[0], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dQdx[1], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dQdy[1], 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(mesh->nx, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(mesh->ny, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(mesh->sJ, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->nu, 0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_gbl(drag, 1, "double", OP_INC),
              op_arg_gbl(lift, 1, "double", OP_INC));

  // Divide by radius of cylinder
  *lift = *lift / 0.05;
  *drag = *drag / 0.05;
}

double Solver::getAvgPressureConvergance() {
  return pressurePoisson->getAverageConvergeIter();
}

double Solver::getAvgViscosityConvergance() {
  return viscosityPoisson->getAverageConvergeIter();
}
