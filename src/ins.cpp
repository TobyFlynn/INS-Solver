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
#include "ins_data.h"
#include "save_solution.h"
#include "poisson.h"
#include "timing.h"
#include "solver.h"

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

int main(int argc, char **argv) {
  op_init(argc, argv, 2);

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

  // Get input from args
  int iter = 1;
  PetscBool found;
  PetscOptionsGetInt(NULL, NULL, "-iter", &iter, &found);

  int save = -1;
  PetscOptionsGetInt(NULL, NULL, "-save", &save, &found);

  int pmethod = 0;
  PetscOptionsGetInt(NULL, NULL, "-pmethod", &pmethod, &found);

  bc_alpha = 0.0;

  Solver *solver = new Solver(filename, pmethod);

  double a0 = 1.0;
  double a1 = 0.0;
  double b0 = 1.0;
  double b1 = 0.0;
  double g0 = 1.0;
  int currentIter = 0;
  double time = 0.0;
  #ifndef INS_MPI
  if(save != -1) {
    save_solution_init("sol.cgns", solver->data);
    export_data_init("data.csv");
  }
  #endif

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
    solver->advection(currentIter % 2, a0, a1, b0, b1, g0, time);
    timer->endAdvection();

    timer->startPressure();
    bool converged = solver->pressure(currentIter % 2, a0, a1, b0, b1, g0, time);
    if(!converged) {
      cout << "******** ERROR ********" << endl;
      cout << "Pressure solve failed to converge, exiting..." << endl;
      cout << "Iteration: " << i << " Time: " << time << endl;
      break;
    }
    timer->endPressure();

    timer->startViscosity();
    converged = solver->viscosity(currentIter % 2, a0, a1, b0, b1, g0, time);
    if(!converged) {
      cout << "******** ERROR ********" << endl;
      cout << "Viscosity solve failed to converge, exiting..." << endl;
      cout << "Iteration: " << i << " Time: " << time << endl;
      break;
    }
    timer->endViscosity();

    currentIter++;
    time += solver->dt;
    #ifndef INS_MPI
    // Calculate drag and lift coefficients + save data
    if(save != -1 && (i + 1) % save == 0) {
      cout << "Iteration: " << i << " Time: " << time << endl;
      timer->startLiftDrag();
      double lift, drag;
      solver->lift_drag_coeff(&lift, &drag, currentIter % 2);
      // export_data("data.csv", i, time, drag, lift, pressurePoisson->getAverageConvergeIter(), viscosityPoisson->getAverageConvergeIter());
      timer->endLiftDrag();

      timer->startSave();
      save_solution_iter("sol.cgns", solver->data, currentIter % 2, (i + 1) / save);
      timer->endSave();
    }
    #endif
  }
  timer->endMainLoop();

  #ifndef INS_MPI
  if(save != -1)
    save_solution_finalise("sol.cgns", solver->data, (iter / save) + 1, solver->dt * save);
  #endif

  // Save solution to CGNS file
  save_solution("end.cgns", solver->data, currentIter % 2);

  timer->endWallTime();
  timer->exportTimings("timings.csv", iter, time);

  cout << "Final time: " << time << endl;
  cout << "Wall time: " << timer->getWallTime() << endl;
  cout << "Solve time: " << timer->getMainLoop() << endl;
  cout << "Time to simulate 1 second: " << timer->getWallTime() / time << endl;
  cout << "Average number of iterations to pressure convergance: " << solver->getAvgPressureConvergance() << endl;
  cout << "Average number of iterations to viscosity convergance: " << solver->getAvgViscosityConvergance() << endl;

  op_timings_to_csv("op2_timings.csv");

  delete solver;
  delete constants;
  delete timer;

  ierr = PetscFinalize();
  // Clean up OP2
  op_exit();
  return ierr;
}
