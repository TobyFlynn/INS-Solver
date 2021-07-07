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

#include "dg_constants.h"
#include "ins_data.h"
#include "save_solution.h"
#include "timing.h"
#include "solver.h"

using namespace std;

Timing *timer;
DGConstants *constants;

extern double reynolds, froude, weber, nu0, nu1, rho0, rho1, ic_u, ic_v;

int main(int argc, char **argv) {
  op_init(argc, argv, 2);

  timer = new Timing();
  timer->startWallTime();
  timer->startSetup();
  constants = new DGConstants();

  char help[] = "Run for i iterations with \"-iter i\"\nSave solution every x iterations with \"-save x\"\n";
  int ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if(ierr) {
    op_printf("Error initialising PETSc\n");
    return ierr;
  }

  // Phi > 0
  nu0 = 1.0;
  rho0 = 1.0;
  // Phi < 0
  nu1 = 1.0;
  rho1 = 1.0;

  ic_u = 0.0;
  ic_v = 0.0;

  double refRho = 1000.0;
  double refMu  = 1e-3;
  double refLen = 0.001;
  double refVel = 1.0;

  // Set Reynolds number
  reynolds = refRho * refVel * refLen / refMu;
  // Set Froude number
  froude = refVel / sqrt(9.8 * refLen);

  // Get input from args
  int iter = 1;
  PetscBool found;
  PetscOptionsGetInt(NULL, NULL, "-iter", &iter, &found);

  int save = -1;
  PetscOptionsGetInt(NULL, NULL, "-save", &save, &found);

  int pre = 0;
  PetscOptionsGetInt(NULL, NULL, "-pre", &pre, &found);

  int problem = 0;
  PetscOptionsGetInt(NULL, NULL, "-problem", &problem, &found);

  int temp = 0;
  PetscOptionsGetInt(NULL, NULL, "-multiphase", &temp, &found);
  bool multiphase = (temp == 1);

  char inputFile[255];
  PetscOptionsGetString(NULL, NULL, "-input", inputFile, 255, &found);
  if(!found) {
    op_printf("Did not specify an input file, use the -input flag\n");
    return -1;
  }
  string filename = string(inputFile);

  char outDir[255];
  string outputDir = "";
  PetscOptionsGetString(NULL, NULL, "-output", outDir, 255, &found);
  if(!found) {
    outputDir = "./";
  } else {
    outputDir = string(outDir);
  }

  if(outputDir.back() != '/') {
    outputDir += "/";
  }

  PetscOptionsGetReal(NULL, NULL, "-mu0", &nu0, &found);
  PetscOptionsGetReal(NULL, NULL, "-mu1", &nu1, &found);
  PetscOptionsGetReal(NULL, NULL, "-rho0", &rho0, &found);
  PetscOptionsGetReal(NULL, NULL, "-rho1", &rho1, &found);

  op_printf("nu0: %g\n", nu0);
  op_printf("nu1: %g\n", nu1);
  op_printf("rho0: %g\n", rho0);
  op_printf("rho1: %g\n", rho1);
  op_printf("reynolds: %g\n", reynolds);
  op_printf("froude: %g\n", froude);

  Solver *solver = new Solver(filename, pre > 0, problem, multiphase);

  double a0 = 1.0;
  double a1 = 0.0;
  double b0 = 1.0;
  double b1 = 0.0;
  double g0 = 1.0;
  int currentIter = 0;
  double time = 0.0;

  if(save != -1) {
    if(multiphase) {
      save_solution_init(outputDir + "sol.cgns", solver->mesh, solver->data, solver->ls);
    } else {
      save_solution_init(outputDir + "sol.cgns", solver->mesh, solver->data, nullptr);
    }
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
    solver->advection(currentIter % 2, a0, a1, b0, b1, g0, time);
    timer->endAdvection();

    timer->startPressure();
    bool converged = solver->pressure(currentIter % 2, a0, a1, b0, b1, g0, time);
    if(!converged) {
      op_printf("******** ERROR ********\n");
      op_printf("Pressure solve failed to converge, exiting...\n");
      op_printf("Iteration: %d Time: %g\n", i, time);
      break;
    }
    timer->endPressure();

    timer->startViscosity();
    converged = solver->viscosity(currentIter % 2, a0, a1, b0, b1, g0, time);
    if(!converged) {
      op_printf("******** ERROR ********\n");
      op_printf("Viscosity solve failed to converge, exiting...\n");
      op_printf("Iteration: %d Time: %g\n", i, time);
      break;
    }
    timer->endViscosity();

    solver->update_surface(currentIter % 2);

    currentIter++;
    time += solver->dt;

    // Save data
    if(save != -1 && (i + 1) % save == 0) {
      op_printf("Iteration: %d Time: %g\n", i, time);

      timer->startSave();
      if(multiphase) {
        save_solution_iter(outputDir + "sol.cgns", solver->mesh, solver->data, currentIter % 2, solver->ls, (i + 1) / save);
      } else {
        save_solution_iter(outputDir + "sol.cgns", solver->mesh, solver->data, currentIter % 2, nullptr, (i + 1) / save);
      }
      timer->endSave();
    }
  }
  timer->endMainLoop();

  if(save != -1)
    save_solution_finalise(outputDir + "sol.cgns", (iter / save) + 1, solver->dt * save);

  // Save solution to CGNS file
  if(multiphase) {
    save_solution(outputDir + "end.cgns", solver->mesh, solver->data, currentIter % 2, solver->ls, time, nu0);
  } else {
    save_solution(outputDir + "end.cgns", solver->mesh, solver->data, currentIter % 2, nullptr, time, nu0);
  }

  timer->endWallTime();
  timer->exportTimings(outputDir + "timings.csv", iter, time);

  op_printf("Final time: %g\n", time);
  op_printf("Wall time: %g\n", timer->getWallTime());
  op_printf("Solve time: %g\n", timer->getMainLoop());
  op_printf("Time to simulate 1 second: %g\n", timer->getWallTime() / time);
  op_printf("Average number of iterations to pressure convergance: %g\n", solver->getAvgPressureConvergance());
  op_printf("Average number of iterations to viscosity convergance: %g\n", solver->getAvgViscosityConvergance());

  string op_out_file = outputDir + "op2_timings.csv";
  op_timings_to_csv(op_out_file.c_str());

  delete solver;
  delete constants;
  delete timer;

  ierr = PetscFinalize();
  // Clean up OP2
  op_exit();
  return ierr;
}
