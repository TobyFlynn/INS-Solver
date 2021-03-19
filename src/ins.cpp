// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <getopt.h>
#include <limits>

#include "constants/all_constants.h"
#include "load_mesh.h"
#include "ins_data.h"
#include "blas_calls.h"
#include "operators.h"
#include "save_solution.h"
#include "poisson.h"

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

#include "kernels/min_max.h"

using namespace std;

// Stuff for parsing command line arguments
extern char *optarg;
extern int  optind, opterr, optopt;
static struct option options[] = {
  {"iter", required_argument, 0, 0},
  {"alpha", required_argument, 0, 0},
  {0,    0,                  0,  0}
};

void advection(INSData *data, int currentInd, double a0, double a1, double b0,
               double b1, double g0, double dt, double t);

void pressure(INSData *data, Poisson *poisson, int currentInd, double a0,
              double a1, double b0, double b1, double g0, double dt, double t);

void viscosity(INSData *data, CubatureData *cubatureData, GaussData *gaussData,
               Poisson *poisson, int currentInd, double a0, double a1, double b0,
               double b1, double g0, double dt, double t);

void get_min_max(INSData *data, double *minQT0, double *minQT1, double *minQTT0,
                 double *minQTT1, double *minQ0, double *minQ1, double *maxQT0,
                 double *maxQT1, double *maxQTT0, double *maxQTT1,
                 double *maxQ0, double *maxQ1, int ind);

int main(int argc, char **argv) {
  string filename = "./cylinder.cgns";
  char help[] = "TODO";
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
  // ic_u = 1.0;
  ic_v = 0.0;

  cout << "gam: " << gam << endl;
  cout << "mu: " << mu << endl;
  cout << "nu: " << nu << endl;

  // Object that holds all sets, maps and dats
  // (along with memory associated with them)
  INSData *data = new INSData();

  auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
    if(x1 == 0.0 && x2 == 0.0) {
      // Inflow
      return 0;
    } else if(x1 == 2.2 && x2 == 2.2) {
      // Outflow
      return 1;
    } else {
      // Wall
      return 2;
    }
  };

  // auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
  //   if(x1 == x2 && x1 < -0.15) {
  //     // Inflow
  //     // cout << "Inflow (" << x1 << "," << y1 << ") and (" << x2 << "," << y2 << ")" << endl;
  //     return 0;
  //   } else if(x1 == x2 && x2 > 1.5) {
  //     // Outflow
  //     // cout << "Outflow (" << x1 << "," << y1 << ") and (" << x2 << "," << y2 << ")" << endl;
  //     return 1;
  //   } else {
  //     // Wall
  //     // cout << "Wall (" << x1 << "," << y1 << ") and (" << x2 << "," << y2 << ")" << endl;
  //     return 2;
  //   }
  // };

  load_mesh(filename.c_str(), data, bcNum);

  // Initialise OP2
  op_init(argc, argv, 2);

  // Initialise all sets, maps and dats
  data->initOP2();

  // Get input from args
  int iter = 1;
  bc_alpha = 0.0;

  int opt_index = 0;
  while(getopt_long_only(argc, argv, "", options, &opt_index) != -1) {
    if(strcmp((char*)options[opt_index].name,"iter") == 0) iter = atoi(optarg);
    if(strcmp((char*)options[opt_index].name,"alpha") == 0) bc_alpha = stod(optarg);
  }

  CubatureData *cubData = new CubatureData(data);
  GaussData *gaussData = new GaussData(data);

  // Set initial conditions
  op_par_loop(set_ic, "set_ic", data->cells,
              op_arg_dat(data->Q[0][0],   -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->Q[0][1],   -1, OP_ID, 15, "double", OP_WRITE));

  cout << "Starting initialisation of pressure Poisson matrix" << endl;
  Poisson *pressurePoisson = new Poisson(data, cubData, gaussData);
  int pressure_dirichlet[] = {1, -1};
  int pressure_neumann[] = {0, 2};
  pressurePoisson->setDirichletBCs(pressure_dirichlet);
  pressurePoisson->setNeumannBCs(pressure_neumann);
  pressurePoisson->createMatrix();
  pressurePoisson->createBCMatrix();
  cout << "Finished initialisation" << endl;
  Poisson *viscosityPoisson = new Poisson(data, cubData, gaussData);
  int viscosity_dirichlet[] = {0, 2};
  int viscosity_neumann[] = {1, -1};
  viscosityPoisson->setDirichletBCs(viscosity_dirichlet);
  viscosityPoisson->setNeumannBCs(viscosity_neumann);
  viscosityPoisson->createMatrix();
  viscosityPoisson->createMassMatrix();
  viscosityPoisson->createBCMatrix();

  double dt = numeric_limits<double>::max();
  op_par_loop(calc_dt, "calc_dt", data->cells,
              op_arg_dat(data->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_MIN));
  // dt = dt * 1e-1;
  // dt = 0.000863006;
  dt = dt / 25.0;
  cout << "dt: " << dt << endl;

  double a0 = 1.0;
  double a1 = 0.0;
  double b0 = 1.0;
  double b1 = 0.0;
  double g0 = 1.0;
  int currentIter = 0;
  double time = 0.0;

  double cpu_1, wall_1, cpu_2, wall_2, cpu_loop_1, wall_loop_1, cpu_loop_2, wall_loop_2;
  double a_time = 0.0;
  double p_time = 0.0;
  double v_time = 0.0;
  op_timers(&cpu_1, &wall_1);

  for(int i = 0; i < iter; i++) {
    if(i == 1) {
      g0 = 1.5;
      a0 = 2.0;
      a1 = -0.5;
      b0 = 2.0;
      b1 = -1.0;
    }
    op_timers(&cpu_loop_1, &wall_loop_1);
    advection(data, currentIter % 2, a0, a1, b0, b1, g0, dt, time);
    op_timers(&cpu_loop_2, &wall_loop_2);
    a_time += wall_loop_2 - wall_loop_1;

    op_timers(&cpu_loop_1, &wall_loop_1);
    pressure(data, pressurePoisson, currentIter % 2, a0, a1, b0, b1, g0, dt, time);
    op_timers(&cpu_loop_2, &wall_loop_2);
    p_time += wall_loop_2 - wall_loop_1;

    op_timers(&cpu_loop_1, &wall_loop_1);
    viscosity(data, cubData, gaussData, viscosityPoisson, currentIter % 2, a0, a1, b0, b1, g0, dt, time);
    op_timers(&cpu_loop_2, &wall_loop_2);
    v_time += wall_loop_2 - wall_loop_1;

    if(i % 100 == 0) {
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

      get_min_max(data, &minQT0, &minQT1, &minQTT0, &minQTT1, &minQ0, &minQ1, &maxQT0,
                  &maxQT1, &maxQTT0, &maxQTT1, &maxQ0, &maxQ1, currentIter % 2);

      cout << "Iter: " << i << endl;
      cout << "QT0: " << minQT0 << " " << maxQT0 << " QT1: " << minQT1 << " " << maxQT1 << endl;
      cout << "QTT0: " << minQTT0 << " " << maxQTT0 << " QTT1: " << minQTT1 << " " << maxQTT1 << endl;
      cout << "Q0: " << minQ0 << " " << maxQ0 << " Q1: " << minQ1 << " " << maxQ1 << endl;
    }
    currentIter++;
    time += dt;
  }
  op_timers(&cpu_2, &wall_2);

  cout << "Final time: " << time << endl;

  cout << "Wall time: " << wall_2 - wall_1 << endl;
  cout << "Time to simulate 1 second: " << (wall_2 - wall_1) / time << endl << endl;

  cout << "Time in advection solve: " << a_time << endl;
  cout << "Time in pressure solve: " << p_time << endl;
  cout << "Time in viscosity solve: " << v_time << endl;

  // Save solution to CGNS file
  double *sol_q0 = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *sol_q1 = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *p_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *pRHS_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *px_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *py_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *utx_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *uty_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *uttx_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *utty_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *visx_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *visy_ptr  = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  op_fetch_data(data->Q[currentIter % 2][0], sol_q0);
  op_fetch_data(data->Q[currentIter % 2][1], sol_q1);
  op_fetch_data(data->dpdx, px_ptr);
  op_fetch_data(data->dpdy, py_ptr);
  op_fetch_data(data->p, p_ptr);
  op_fetch_data(data->pRHS, pRHS_ptr);
  op_fetch_data(data->QT[0], utx_ptr);
  op_fetch_data(data->QT[1], uty_ptr);
  op_fetch_data(data->QTT[0], uttx_ptr);
  op_fetch_data(data->QTT[1], utty_ptr);
  op_fetch_data(data->visRHS[0], visx_ptr);
  op_fetch_data(data->visRHS[1], visy_ptr);
  // save_solution_cell("cylinder.cgns", op_get_size(data->nodes), op_get_size(data->cells),
  //               sol_q0, p_ptr, data->cgnsCells);

  // op_fetch_data(data->p, sol_q0);
  // op_fetch_data(data->Q[currentIter % 2][1], sol_q1);
  // save_solution("cylinder.cgns", op_get_size(data->nodes), op_get_size(data->cells),
  //               sol_q0, sol_q1, p_ptr, data->cgnsCells);

  save_solution_t("cylinder.cgns", op_get_size(data->nodes), op_get_size(data->cells),
                  sol_q0, sol_q1, p_ptr, pRHS_ptr, px_ptr, py_ptr, utx_ptr, uty_ptr, uttx_ptr, utty_ptr, visx_ptr, visy_ptr, data->cgnsCells);

  free(sol_q0);
  free(sol_q1);
  free(p_ptr);
  free(pRHS_ptr);
  free(px_ptr);
  free(py_ptr);
  free(utx_ptr);
  free(uty_ptr);
  free(uttx_ptr);
  free(utty_ptr);
  free(visx_ptr);
  free(visy_ptr);

  // op_fetch_data_hdf5_file(data->Q[currentIter % 2][0], "sol.h5");
  // op_fetch_data_hdf5_file(data->Q[currentIter % 2][1], "sol.h5");

  // Clean up OP2
  op_exit();

  delete pressurePoisson;
  delete viscosityPoisson;
  delete gaussData;
  delete cubData;
  delete data;

  ierr = PetscFinalize();
  return ierr;
}

void advection(INSData *data, int currentInd, double a0, double a1, double b0,
               double b1, double g0, double dt, double t) {
  op_par_loop(advection_flux, "advection_flux", data->cells,
              op_arg_dat(data->Q[currentInd][0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->F[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->F[3], -1, OP_ID, 15, "double", OP_WRITE));

  div(data, data->F[0], data->F[1], data->N[currentInd][0]);
  div(data, data->F[2], data->F[3], data->N[currentInd][1]);

  op_par_loop(advection_faces, "advection_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][0], -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->Q[currentInd][1], -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->exQ[0], -2, data->edge2cells, 15, "double", OP_INC),
              op_arg_dat(data->exQ[1], -2, data->edge2cells, 15, "double", OP_INC));

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

void pressure(INSData *data, Poisson *poisson, int currentInd, double a0, double a1, double b0,
              double b1, double g0, double dt, double t) {
  div(data, data->QT[0], data->QT[1], data->divVelT);
  curl(data, data->Q[currentInd][0], data->Q[currentInd][1], data->curlVel);
  grad(data, data->curlVel, data->gradCurlVel[0], data->gradCurlVel[1]);

  // Apply boundary conditions
  // May need to change in future if non-constant boundary conditions used
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

  // op_fetch_data_hdf5_file(data->dPdN[currentInd], "dpdn.h5");
  // op_fetch_data_hdf5_file(data->divVelT, "div.h5");

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

  // op_fetch_data_hdf5_file(data->dPdN[(currentInd + 1) % 2], "dpdn2.h5");

  pressure_rhs_blas(data, currentInd);

  poisson->setBCValues(data->dirichletBC);
  poisson->solve(data->pRHS, data->p);

  grad(data, data->p, data->dpdx, data->dpdy);

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
}

void viscosity(INSData *data, CubatureData *cubatureData, GaussData *gaussData,
               Poisson *poisson, int currentInd, double a0, double a1, double b0,
               double b1, double g0, double dt, double t) {
  double time = t + dt;

  op_par_loop(viscosity_bc, "viscosity_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&time, 1, "double", OP_READ),
              op_arg_dat(gaussData->x, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gaussData->y, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(data->visBC[0], 0, data->bedge2cells, 21, "double", OP_INC),
              op_arg_dat(data->visBC[1], 0, data->bedge2cells, 21, "double", OP_INC));

  viscosity_rhs_blas(data, cubatureData);

  double factor = g0 / (nu * dt);
  op_par_loop(viscosity_rhs, "viscosity_rhs", data->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->visRHS[0], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->visRHS[1], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->visBC[0], -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(data->visBC[1], -1, OP_ID, 21, "double", OP_RW));

  poisson->setBCValues(data->visBC[0]);
  poisson->solve(data->visRHS[0], data->Q[(currentInd + 1) % 2][0], true, factor);

  poisson->setBCValues(data->visBC[1]);
  poisson->solve(data->visRHS[1], data->Q[(currentInd + 1) % 2][1], true, factor);

  op_par_loop(viscosity_reset_bc, "viscosity_reset_bc", data->cells,
              op_arg_dat(data->visBC[0], -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(data->visBC[1], -1, OP_ID, 21, "double", OP_WRITE));
}

void get_min_max(INSData *data, double *minQT0, double *minQT1, double *minQTT0,
                 double *minQTT1, double *minQ0, double *minQ1, double *maxQT0,
                 double *maxQT1, double *maxQTT0, double *maxQTT1,
                 double *maxQ0, double *maxQ1, int ind) {
  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(minQT0, 1, "double", OP_MIN),
              op_arg_gbl(maxQT0, 1, "double", OP_MAX),
              op_arg_dat(data->QT[0], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(minQT1, 1, "double", OP_MIN),
              op_arg_gbl(maxQT1, 1, "double", OP_MAX),
              op_arg_dat(data->QT[1], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(minQTT0, 1, "double", OP_MIN),
              op_arg_gbl(maxQTT0, 1, "double", OP_MAX),
              op_arg_dat(data->QTT[0], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(minQTT1, 1, "double", OP_MIN),
              op_arg_gbl(maxQTT1, 1, "double", OP_MAX),
              op_arg_dat(data->QTT[1], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(minQ0, 1, "double", OP_MIN),
              op_arg_gbl(maxQ0, 1, "double", OP_MAX),
              op_arg_dat(data->Q[(ind + 1) % 2][0], -1, OP_ID, 15, "double", OP_READ));

  op_par_loop(min_max, "min_max", data->cells,
              op_arg_gbl(minQ1, 1, "double", OP_MIN),
              op_arg_gbl(maxQ1, 1, "double", OP_MAX),
              op_arg_dat(data->Q[(ind + 1) % 2][1], -1, OP_ID, 15, "double", OP_READ));
}
