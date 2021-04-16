// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

#include "constants/all_constants.h"
#include "constants.h"
#include "load_mesh.h"
#include "ins_data.h"
#include "blas_calls.h"
#include "operators.h"
#include "save_solution.h"
#include "poisson.h"
#include "timing.h"

// Kernels
#include "kernels/set_ic.h"
#include "kernels/calc_dt.h"

#include "kernels/advection_flux.h"
#include "kernels/advection_faces.h"
#include "kernels/advection_bc.h"
#include "kernels/advection_numerical_flux.h"
#include "kernels/advection_intermediate_vel.h"

#include "kernels/pressure_bc.h"
#include "kernels/pressure_rhs.h"
#include "kernels/pressure_update_vel.h"

#include "kernels/viscosity_rhs.h"
#include "kernels/viscosity_reset_bc.h"
#include "kernels/viscosity_bc.h"

#include "kernels/lift_drag.h"
#include "kernels/min_max.h"

using namespace std;

void export_data_init(string filename) {
  ofstream file(filename);

  // Create first row of csv file (headers of columns)
  file << "Iteration" << ",";
  file << "Time" << ",";
  file << "Drag Coefficient" << ",";
  file << "Lift Coefficient" << ",";
  file << "Avg. Pressure Convergance" << ",";
  file << "Avg. Viscosity Convergance" << endl;

  file.close();
}

void export_data(string filename, int iter, double time, double drag,
                 double lift, double avgPr, double avgVis) {
  ofstream file(filename, ios::app);

  // Create first row of csv file (headers of columns)
  file << to_string(iter) << ",";
  file << to_string(time) << ",";
  file << to_string(drag) << ",";
  file << to_string(lift) << ",";
  file << to_string(avgPr) << ",";
  file << to_string(avgVis) << endl;

  file.close();
}

Timing *timer;
Constants *constants;

void advection(INSData *data, int currentInd, double a0, double a1, double b0,
               double b1, double g0, double dt, double t);

bool pressure(INSData *data, Poisson *poisson, int currentInd, double a0,
              double a1, double b0, double b1, double g0, double dt, double t);

bool viscosity(INSData *data, CubatureData *cubatureData, GaussData *gaussData,
               Poisson *poisson, int currentInd, double a0, double a1, double b0,
               double b1, double g0, double dt, double t);

void lift_drag_coeff(INSData *data, double *lift, double *drag, int ind);

void print_min_max(INSData *data, int ind);

int main(int argc, char **argv) {
  timer = new Timing();
  timer->startWallTime();
  timer->startSetup();
  constants = new Constants();

  string filename = "./cylinder.cgns";
  char help[] = "Run for i iterations with \"-iter i\"\nSave solution every x iterations with \"-save x\"\n";
  int ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if(ierr) {
    cout << "Error initialising PETSc" << endl;
    return ierr;
  }

  gam = 1.4;
  mu = 1e-2;
  nu = 1e-3;
  bc_u = 1e-6;
  bc_v = 0.0;
  ic_u = 0.0;
  ic_v = 0.0;

  cout << "gam: " << gam << endl;
  cout << "mu: " << mu << endl;
  cout << "nu: " << nu << endl;

  // Object that holds all sets, maps and dats (along with memory associated with them)
  INSData *data = new INSData();

  // Lamda used to identify the type of boundary edges
  auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
    if(x1 == 0.0 && x2 == 0.0) {
      // Inflow
      return 0;
    } else if(x1 == 2.2 && x2 == 2.2) {
      // Outflow
      return 1;
    } else if(x1 > 0.1 && x2 > 0.1 && x1 < 0.3 && x2 < 0.3
              && y1 > 0.1 && y2 > 0.1 && y1 < 0.3 && y2 < 0.3) {
      // Cylinder Wall
      return 2;
    } else {
      // Top/Bottom Wall
      return 3;
    }
  };

  load_mesh(filename.c_str(), data, bcNum);

  // Initialise OP2
  op_init(argc, argv, 2);

  // Initialise all sets, maps and dats
  data->initOP2();

  // Get input from args
  int iter = 1;
  PetscBool found;
  PetscOptionsGetInt(NULL, NULL, "-iter", &iter, &found);

  int save = -1;
  PetscOptionsGetInt(NULL, NULL, "-save", &save, &found);

  bc_alpha = 0.0;

  CubatureData *cubData = new CubatureData(data);
  GaussData *gaussData = new GaussData(data);

  // Set initial conditions
  op_par_loop(set_ic, "set_ic", data->cells,
              op_arg_dat(data->Q[0][0],   -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->Q[0][1],   -1, OP_ID, 15, "double", OP_WRITE));

  // Initialise Poisson solvers
  // Poisson_M *pressurePoisson = new Poisson_M(data, cubData, gaussData);
  // int pressure_dirichlet[] = {1, -1, -1};
  // int pressure_neumann[] = {0, 2, 3};
  // pressurePoisson->setDirichletBCs(pressure_dirichlet);
  // pressurePoisson->setNeumannBCs(pressure_neumann);
  // pressurePoisson->createMatrix();
  // pressurePoisson->createBCMatrix();
  Poisson_MF *pressurePoisson = new Poisson_MF(data, cubData, gaussData);
  int pressure_dirichlet[] = {1, -1, -1};
  int pressure_neumann[] = {0, 2, 3};
  pressurePoisson->setDirichletBCs(pressure_dirichlet);
  pressurePoisson->setNeumannBCs(pressure_neumann);
  pressurePoisson->createBCMatrix();
  // Poisson_M *viscosityPoisson = new Poisson_M(data, cubData, gaussData);
  // int viscosity_dirichlet[] = {0, 2, 3};
  // int viscosity_neumann[] = {1, -1, -1};
  // viscosityPoisson->setDirichletBCs(viscosity_dirichlet);
  // viscosityPoisson->setNeumannBCs(viscosity_neumann);
  // viscosityPoisson->createMatrix();
  // viscosityPoisson->createMassMatrix();
  // viscosityPoisson->createBCMatrix();
  Poisson_MF *viscosityPoisson = new Poisson_MF(data, cubData, gaussData);
  int viscosity_dirichlet[] = {0, 2, 3};
  int viscosity_neumann[] = {1, -1, -1};
  viscosityPoisson->setDirichletBCs(viscosity_dirichlet);
  viscosityPoisson->setNeumannBCs(viscosity_neumann);
  viscosityPoisson->createBCMatrix();

  double dt = numeric_limits<double>::max();
  op_par_loop(calc_dt, "calc_dt", data->cells,
              op_arg_dat(data->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_MIN));
  dt = dt / 25.0;
  cout << "dt: " << dt << endl;

  double a0 = 1.0;
  double a1 = 0.0;
  double b0 = 1.0;
  double b1 = 0.0;
  double g0 = 1.0;
  int currentIter = 0;
  double time = 0.0;

  if(save != -1) {
    save_solution_init("sol.cgns", data);
    export_data_init("data.csv");
  }

  timer->endSetup();
  timer->startMainLoop();

  for(int i = 0; i < iter; i++) {
    // Switch from forwards Euler time integration to second-order Adams-Bashford after first iteration
    if(i == 1) {
      g0 = 1.5;
      a0 = 2.0;
      a1 = -0.5;
      b0 = 2.0;
      b1 = -1.0;
    }
    timer->startAdvection();
    advection(data, currentIter % 2, a0, a1, b0, b1, g0, dt, time);
    timer->endAdvection();

    timer->startPressure();
    bool converged = pressure(data, pressurePoisson, currentIter % 2, a0, a1, b0, b1, g0, dt, time);
    if(!converged) {
      cout << "******** ERROR ********" << endl;
      cout << "Pressure solve failed to converge, exiting..." << endl;
      cout << "Iteration: " << i << " Time: " << time << endl;
      break;
    }
    timer->endPressure();

    timer->startViscosity();
    converged = viscosity(data, cubData, gaussData, viscosityPoisson, currentIter % 2, a0, a1, b0, b1, g0, dt, time);
    if(!converged) {
      cout << "******** ERROR ********" << endl;
      cout << "Viscosity solve failed to converge, exiting..." << endl;
      cout << "Iteration: " << i << " Time: " << time << endl;
      break;
    }
    timer->endViscosity();

    currentIter++;
    time += dt;

    // Calculate drag and lift coefficients + save data
    if(save != -1 && (i + 1) % save == 0) {
      // print_min_max(data, currentIter % 2);
      cout << "Iteration: " << i << " Time: " << time << endl;
      timer->startLiftDrag();
      double lift, drag;
      lift_drag_coeff(data, &lift, &drag, currentIter % 2);
      export_data("data.csv", i, time, drag, lift, pressurePoisson->getAverageConvergeIter(), viscosityPoisson->getAverageConvergeIter());
      timer->endLiftDrag();

      timer->startSave();
      save_solution_iter("sol.cgns", data, currentIter % 2, (i + 1) / save);
      timer->endSave();
    }
  }
  timer->endMainLoop();

  if(save != -1)
    save_solution_finalise("sol.cgns", data, (iter / save) + 1, dt * save);

  // Save solution to CGNS file
  save_solution("end.cgns", data, currentIter % 2);

  timer->endWallTime();
  timer->exportTimings("timings.csv", iter, time);

  cout << "Final time: " << time << endl;
  cout << "Wall time: " << timer->getWallTime() << endl;
  cout << "Time to simulate 1 second: " << timer->getWallTime() / time << endl;

  op_timings_to_csv("op2_timings.csv");

  // Clean up OP2
  op_exit();

  delete pressurePoisson;
  delete viscosityPoisson;
  delete gaussData;
  delete cubData;
  delete data;
  delete constants;
  delete timer;

  ierr = PetscFinalize();
  return ierr;
}

void advection(INSData *data, int currentInd, double a0, double a1, double b0,
               double b1, double g0, double dt, double t) {
  // Calculate flux values
  op_par_loop(advection_flux, "advection_flux", data->cells,
              op_arg_dat(data->Q[currentInd][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->F[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[3], -1, OP_ID, 15, "double", OP_WRITE));

  div(data, data->F[0], data->F[1], data->N[currentInd][0]);
  div(data, data->F[2], data->F[3], data->N[currentInd][1]);

  // Exchange values on edges between elements
  op_par_loop(advection_faces, "advection_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][0], -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->exQ[0], -2, data->edge2cells, 15, "double", OP_INC),
              op_arg_dat(data->exQ[1], -2, data->edge2cells, 15, "double", OP_INC));

  // Enforce BCs
  op_par_loop(advection_bc, "advection_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_dat(data->x, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->y, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][0], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->exQ[0], 0, data->bedge2cells, 15, "double", OP_INC),
              op_arg_dat(data->exQ[1], 0, data->bedge2cells, 15, "double", OP_INC));

  // Calculate numberical flux across edges
  op_par_loop(advection_numerical_flux, "advection_numerical_flux", data->cells,
              op_arg_dat(data->fscale,  -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->nx,      -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ny,      -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->exQ[0],  -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->exQ[1],  -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->flux[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->flux[1], -1, OP_ID, 15, "double", OP_WRITE));

  advection_lift_blas(data, currentInd);

  // Calculate the intermediate velocity values
  op_par_loop(advection_intermediate_vel, "advection_intermediate_vel", data->cells,
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

bool pressure(INSData *data, Poisson *poisson, int currentInd, double a0, double a1, double b0,
              double b1, double g0, double dt, double t) {
  timer->startPressureSetup();
  div(data, data->QT[0], data->QT[1], data->divVelT);
  curl(data, data->Q[currentInd][0], data->Q[currentInd][1], data->curlVel);
  grad(data, data->curlVel, data->gradCurlVel[0], data->gradCurlVel[1]);

  // Apply pressure boundary conditions
  op_par_loop(pressure_bc, "pressure_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_dat(data->x, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->y, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->nx, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->ny, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->N[currentInd][0], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->N[currentInd][1], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->gradCurlVel[0], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->gradCurlVel[1], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], 0, data->bedge2cells, 15, "double", OP_INC));

  // Calculate RHS of pressure solve
  op_par_loop(pressure_rhs, "pressure_rhs", data->cells,
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&g0, 1, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sJ, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->divVelT, -1, OP_ID, 15, "double", OP_RW));

  pressure_rhs_blas(data, currentInd);
  timer->endPressureSetup();

  // Call PETSc linear solver
  timer->startPressureLinearSolve();
  poisson->setBCValues(data->zeroBC);
  bool converged = poisson->solve(data->pRHS, data->p);
  timer->endPressureLinearSolve();

  // Calculate gradient of pressure
  grad(data, data->p, data->dpdx, data->dpdy);

  // Calculate new velocity intermediate values
  double factor = dt / g0;
  op_par_loop(pressure_update_vel, "pressure_update_vel", data->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(data->dpdx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->dpdy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->QT[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->QT[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->QTT[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->QTT[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, 15, "double", OP_WRITE));

  return converged;
}

bool viscosity(INSData *data, CubatureData *cubatureData, GaussData *gaussData,
               Poisson *poisson, int currentInd, double a0, double a1, double b0,
               double b1, double g0, double dt, double t) {
  timer->startViscositySetup();
  double time = t + dt;
  // Get BCs for viscosity solve
  op_par_loop(viscosity_bc, "viscosity_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&time, 1, "double", OP_READ),
              op_arg_dat(gaussData->x, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gaussData->y, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(data->visBC[0], 0, data->bedge2cells, 21, "double", OP_INC),
              op_arg_dat(data->visBC[1], 0, data->bedge2cells, 21, "double", OP_INC));

  // Set up RHS for viscosity solve
  viscosity_rhs_blas(data, cubatureData);

  double factor = g0 / (nu * dt);
  op_par_loop(viscosity_rhs, "viscosity_rhs", data->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->visRHS[0], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->visRHS[1], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->visBC[0], -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(data->visBC[1], -1, OP_ID, 21, "double", OP_RW));

  timer->endViscositySetup();

  // Call PETSc linear solver
  timer->startViscosityLinearSolve();
  poisson->setBCValues(data->visBC[0]);
  bool convergedX = poisson->solve(data->visRHS[0], data->Q[(currentInd + 1) % 2][0], true, factor);

  poisson->setBCValues(data->visBC[1]);
  bool convergedY = poisson->solve(data->visRHS[1], data->Q[(currentInd + 1) % 2][1], true, factor);
  timer->endViscosityLinearSolve();

  // Reset BC dats ready for next iteration
  op_par_loop(viscosity_reset_bc, "viscosity_reset_bc", data->cells,
              op_arg_dat(data->visBC[0], -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(data->visBC[1], -1, OP_ID, 21, "double", OP_WRITE));

  return convergedX && convergedY;
}

// Function to calculate lift and drag coefficients of the cylinder
void lift_drag_coeff(INSData *data, double *lift, double *drag, int ind) {
  *lift = 0.0;
  *drag = 0.0;

  grad(data, data->Q[(ind + 1) % 2][0], data->dQdx[0], data->dQdy[0]);
  grad(data, data->Q[(ind + 1) % 2][1], data->dQdx[1], data->dQdy[1]);

  op_par_loop(lift_drag, "lift_drag", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->p, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dQdx[0], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dQdy[0], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dQdx[1], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dQdy[1], 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->nx, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->ny, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->sJ, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_gbl(drag, 1, "double", OP_INC),
              op_arg_gbl(lift, 1, "double", OP_INC));

  // Divide by radius of cylinder
  *lift = *lift / 0.05;
  *drag = *drag / 0.05;
}

// Function for printing the min and max values of various OP2 dats
void print_min_max(INSData *data, int ind) {
  double minQT0 = numeric_limits<double>::max();
  double minQT1 = numeric_limits<double>::max();
  double minQTT0 = numeric_limits<double>::max();
  double minQTT1 = numeric_limits<double>::max();
  double minQ0 = numeric_limits<double>::max();
  double minQ1 = numeric_limits<double>::max();
  double maxQT0 = 0.0;
  double maxQT1 = 0.0;
  double maxQTT0 = 0.0;
  double maxQTT1 = 0.0;
  double maxQ0 = 0.0;
  double maxQ1 = 0.0;

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(&minQT0, 1, "double", OP_MIN),
              op_arg_gbl(&maxQT0, 1, "double", OP_MAX),
              op_arg_dat(data->QT[0], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(&minQT1, 1, "double", OP_MIN),
              op_arg_gbl(&maxQT1, 1, "double", OP_MAX),
              op_arg_dat(data->QT[1], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(&minQTT0, 1, "double", OP_MIN),
              op_arg_gbl(&maxQTT0, 1, "double", OP_MAX),
              op_arg_dat(data->QTT[0], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(&minQTT1, 1, "double", OP_MIN),
              op_arg_gbl(&maxQTT1, 1, "double", OP_MAX),
              op_arg_dat(data->QTT[1], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(&minQ0, 1, "double", OP_MIN),
              op_arg_gbl(&maxQ0, 1, "double", OP_MAX),
              op_arg_dat(data->Q[ind][0], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(&minQ1, 1, "double", OP_MIN),
              op_arg_gbl(&maxQ1, 1, "double", OP_MAX),
              op_arg_dat(data->Q[ind][1], -1, OP_ID, 15, "double", OP_READ));

  cout << "QT0: " << minQT0 << " " << maxQT0 << " QT1: " << minQT1 << " " << maxQT1 << endl;
  cout << "QTT0: " << minQTT0 << " " << maxQTT0 << " QTT1: " << minQTT1 << " " << maxQTT1 << endl;
  cout << "Q0: " << minQ0 << " " << maxQ0 << " Q1: " << minQ1 << " " << maxQ1 << endl;
}
