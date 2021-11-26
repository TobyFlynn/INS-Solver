#include "solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_constants.h"
#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_operators.h"
#include "load_mesh.h"
#include "timing.h"

extern Timing *timer;
extern double dt, reynolds, refVel;

using namespace std;

Solver::Solver(std::string filename, int prob) {
  problem = prob;

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
  ls = new LS(mesh, data);
  pressurePoisson = new PressureSolve(mesh, data);
  viscosityPoisson = new ViscositySolve(mesh, data);

  pressurePoisson->setDirichletBCs(pressure_dirichlet);
  pressurePoisson->setNeumannBCs(pressure_neumann);
  viscosityPoisson->setDirichletBCs(viscosity_dirichlet);
  viscosityPoisson->setNeumannBCs(viscosity_neumann);

  op_partition("PARMETIS", "KWAY", mesh->cells, mesh->edge2cells, NULL);

  mesh->init();
  data->init();
  ls->init();
  pressurePoisson->init();
  viscosityPoisson->init();

  // Set initial conditions
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[0][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[0][1], -1, OP_ID, DG_NP, "double", OP_WRITE));

  dt = numeric_limits<double>::max();
  op_par_loop(calc_dt, "calc_dt", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_MIN));
  dt = dt / (DG_ORDER * DG_ORDER * refVel);
  op_printf("dt: %g\n", dt);
}

Solver::~Solver() {
  delete viscosityPoisson;
  delete pressurePoisson;
  delete ls;
  delete data;
  delete mesh;
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
  op_par_loop(advection_bc, "advection_bc", mesh->bedges,
              op_arg_dat(mesh->order,       0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->x,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->y,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[0],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[1],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->flux[0],   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->flux[1],   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  op2_gemv(mesh, true, 1.0, DGConstants::GAUSS_INTERP, data->flux[0], 1.0, data->N[currentInd][0]);
  op2_gemv(mesh, true, 1.0, DGConstants::GAUSS_INTERP, data->flux[1], 1.0, data->N[currentInd][1]);
  // op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[0], 1.0, data->N[currentInd][0]);
  // op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[1], 1.0, data->N[currentInd][1]);

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
  timer->startPressureSetup();
  op_par_loop(save_order, "save_order", mesh->cells,
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);
  dats_to_update.push_back(data->QT[0]);
  dats_to_update.push_back(data->QT[1]);
  dats_to_update.push_back(data->N[0][0]);
  dats_to_update.push_back(data->N[0][1]);
  dats_to_update.push_back(data->N[1][0]);
  dats_to_update.push_back(data->N[1][1]);
  mesh->update_order(data->new_order, dats_to_update);

  div(mesh, data->QT[0], data->QT[1], data->divVelT);
  curl(mesh, data->Q[currentInd][0], data->Q[currentInd][1], data->curlVel);
  grad(mesh, data->curlVel, data->gradCurlVel[0], data->gradCurlVel[1]);

  // Apply pressure boundary conditions
  op_par_loop(pressure_bc, "pressure_bc", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->order, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->nx, 0, mesh->bedge2cells, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(mesh->ny, 0, mesh->bedge2cells, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(data->N[currentInd][0], 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->N[currentInd][1], 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->gradCurlVel[0], 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->gradCurlVel[1], 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], 0, mesh->bedge2cells, 3 * DG_NPF, "double", OP_INC));

  if(problem == 1) {
    // TODO - update for p-adaptivity
    op_par_loop(pressure_bc2, "pressure_bc2", mesh->bedges,
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_gbl(&t, 1, "double", OP_READ),
                op_arg_gbl(&problem, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->prBC, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));
  }

  // Calculate RHS of pressure solve
  op_par_loop(pressure_rhs, "pressure_rhs", mesh->cells,
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], -1, OP_ID, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, 3 * DG_NPF, "double", OP_RW),
              op_arg_dat(data->divVelT, -1, OP_ID, DG_NP, "double", OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, data->dPdN[(currentInd + 1) % 2], 1.0, data->divVelT);
  op2_gemv(mesh, false, 1.0, DGConstants::MASS, data->divVelT, 0.0, data->pRHS);
  timer->endPressureSetup();

  // Call PETSc linear solver
  timer->startPressureLinearSolve();
  pressurePoisson->setup();
  pressurePoisson->setBCValues(data->prBC);
  bool converged = pressurePoisson->solve(data->pRHS, data->p);
  timer->endPressureLinearSolve();

  // Calculate gradient of pressure
  grad(mesh, data->p, data->dpdx, data->dpdy);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, data->p, 0.0, data->gP);

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(data->pFluxX, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(data->pFluxY, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  op_par_loop(pressure_grad_flux, "pressure_grad_flux", mesh->edges,
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gP,        -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->pFluxX,    -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->pFluxY,    -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  // op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, data->pFluxX, 1.0, data->dpdx);
  // op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, data->pFluxY, 1.0, data->dpdy);
  op2_gemv(mesh, false, -1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->pFluxX, 1.0, data->dpdx);
  op2_gemv(mesh, false, -1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->pFluxY, 1.0, data->dpdy);

  // Calculate new velocity intermediate values
  op_par_loop(pressure_update_vel, "pressure_update_vel", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(data->dpdx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->dpdy, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QT[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QT[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->QTT[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, 3 * DG_NPF, "double", OP_WRITE),
              op_arg_dat(data->prBC, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  return converged;
}

bool Solver::viscosity(int currentInd, double a0, double a1, double b0,
                       double b1, double g0, double t) {
  timer->startViscositySetup();
  double time = t + dt;

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(data->visBC[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(data->visBC[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Get BCs for viscosity solve
  op_par_loop(viscosity_bc, "viscosity_bc", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&time, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->visBC[0], 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->visBC[1], 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  // Set up RHS for viscosity solve
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

  double factor = reynolds / dt;
  op_par_loop(viscosity_rhs, "viscosity_rhs", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(data->visRHS[0], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(data->visRHS[1], -1, OP_ID, DG_NP, "double", OP_RW));

  timer->endViscositySetup();

  // Call PETSc linear solver
  timer->startViscosityLinearSolve();
  factor = g0 * reynolds / dt;
  viscosityPoisson->setup(factor);
  viscosityPoisson->setBCValues(data->visBC[0]);
  bool convergedX = viscosityPoisson->solve(data->visRHS[0], data->Q[(currentInd + 1) % 2][0]);

  viscosityPoisson->setBCValues(data->visBC[1]);
  bool convergedY = viscosityPoisson->solve(data->visRHS[1], data->Q[(currentInd + 1) % 2][1]);
  timer->endViscosityLinearSolve();

  op_par_loop(ls_update_order, "ls_update_order", mesh->cells,
              op_arg_dat(mesh->order,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&ls->alpha,       1, "double", OP_READ),
              op_arg_dat(ls->s,           -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);
  dats_to_update.push_back(ls->s);

  mesh->update_order(data->new_order, dats_to_update);

  return convergedX && convergedY;
}

void Solver::update_surface(int currentInd) {
  timer->startSurface();
  ls->setVelField(data->Q[(currentInd + 1) % 2][0], data->Q[(currentInd + 1) % 2][1]);
  ls->step(dt);
  timer->endSurface();
}

double Solver::getAvgPressureConvergance() {
  return pressurePoisson->getAverageConvergeIter();
}

double Solver::getAvgViscosityConvergance() {
  return viscosityPoisson->getAverageConvergeIter();
}
