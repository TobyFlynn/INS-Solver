#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>
#include <set>
#include <sstream>

#include "petscvec.h"
#include "petscksp.h"

#include "dg_global_constants/dg_global_constants_2d.h"
#include "dg_abort.h"

#include "timing.h"
#include "config.h"
#include "solvers/2d/ins_solver.h"
#include "solvers/2d/ins_temperature_solver.h"
#include "solvers/2d/mp_ins_solver.h"
#include "drivers/2d/ls_driver.h"
#include "drivers/2d/compressible_euler_driver.h"
#include "measurements/2d/lift_drag_cylinder.h"
#include "measurements/2d/l2_error_vortex.h"
#include "measurements/2d/min_max_interface.h"
#include "measurements/2d/min_max_pressure.h"
#include "measurements/2d/ls_advec_error.h"
#include "measurements/2d/mass_of_phases.h"
#include "measurements/2d/l2_euler_vortex.h"

Timing *timer;
Config *config;

using namespace std;

// Global constants
DG_FP r_ynolds, mu0, mu1, rho0, rho1, gamma_e, weber, froude, coeff_thermal_expan, peclet;

// Set up measurements
void add_measurements(SimulationDriver *driver, vector<Measurement2D*> &measurements, const DG_FP refMu) {
  // Get list of measurements to take
  // Options are: lift_drag, l2_vortex, min_max_interface
  string mes_tmp = "none";
  config->getStr("io", "measurements", mes_tmp);
  if(mes_tmp != "none") {
    set<string> measurements_to_take;
    stringstream tmp_ss(mes_tmp);
    string val_str;
    while(getline(tmp_ss, val_str, ',')) {
      measurements_to_take.insert(val_str);
    }

    for(auto &measurement : measurements_to_take) {
      if(measurement == "lift_drag") {
        LiftDragCylinder2D *lift_drag = new LiftDragCylinder2D(driver, 7.0, 7.0, 9.0, 9.0);
        measurements.push_back(lift_drag);
      } else if(measurement == "l2_vortex") {
        L2ErrorVortex2D *l2_error = new L2ErrorVortex2D(driver);
        measurements.push_back(l2_error);
      } else if(measurement == "min_max_interface") {
        MinMaxInterface2D *min_max = new MinMaxInterface2D(driver);
        measurements.push_back(min_max);
      } else if(measurement == "min_max_pressure") {
        MinMaxPressure2D *min_max_pr = new MinMaxPressure2D(driver, 1);
        measurements.push_back(min_max_pr);
      } else if(measurement == "l2_ls_advec") {
        LSAdvecError *l2_error = new LSAdvecError(driver, 1);
        measurements.push_back(l2_error);
      } else if(measurement == "mass_of_phases") {
        MassOfPhases2D *mass_error = new MassOfPhases2D(driver);
        measurements.push_back(mass_error);
      } else if(measurement == "l2_euler_vortex") {
        L2EulerVortex2D *l2_error = new L2EulerVortex2D(driver);
        measurements.push_back(l2_error);
      } else {
        dg_abort("Unrecognised measurement: " + measurement);
      }
    }
  }
}

int main(int argc, char **argv) {
  op_init(argc, argv, 1);

  timer = new Timing();

  char help[] = "Run for i iterations with \"-iter i\"\nSave solution every x iterations with \"-save x\"\n";
  int ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if(ierr) {
    op_printf("Error initialising PETSc\n");
    return ierr;
  }

  PetscBool found;
  char inputFile[255];
  PetscOptionsGetString(NULL, NULL, "-input", inputFile, 255, &found);
  if(!found) {
    op_printf("Did not specify an input file, use the -input flag\n");
    return -1;
  }
  string filename = string(inputFile);

  char configFile[255];
  PetscOptionsGetString(NULL, NULL, "-config", configFile, 255, &found);
  if(!found) {
    op_printf("Did not specify a configuration file, use the -config flag\n");
    return -1;
  }
  config = new Config(string(configFile));

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

  bool resuming = false;
  int resumeIter = 0;
  char checkFile[255];
  string checkpointFile;
  PetscOptionsGetString(NULL, NULL, "-checkpoint", checkFile, 255, &found);
  if(found) {
    checkpointFile = string(checkFile);
    resuming = true;
    op_get_const_hdf5("iter", 1, "int", (char *)&resumeIter, checkpointFile.c_str());
  }

  int iter = 1;
  config->getInt("simulation-constants", "iter", iter);

  int save = -1;
  config->getInt("io", "save", save);
  int checkpoint = -1;
  config->getInt("io", "checkpoint", checkpoint);

  mu0  = 1.0;
  mu1  = 1.0;
  rho0 = 1.0;
  rho1 = 1.0;
  config->getDouble("fluid-constants", "mu0", mu0);
  config->getDouble("fluid-constants", "mu1", mu1);
  config->getDouble("fluid-constants", "rho0", rho0);
  config->getDouble("fluid-constants", "rho1", rho1);

  DG_FP refRho = 1.0;
  DG_FP refVel = 1.0;
  DG_FP refLen = 0.005;
  DG_FP refMu  = 1.0e-5;
  DG_FP refSurfTen = 70.0 * 1e-3;
  DG_FP refGrav = 9.8;
  config->getDouble("fluid-constants", "refRho", refRho);
  config->getDouble("fluid-constants", "refVel", refVel);
  config->getDouble("fluid-constants", "refLen", refLen);
  config->getDouble("fluid-constants", "refMu", refMu);
  config->getDouble("fluid-constants", "refSurfTen", refSurfTen);
  config->getDouble("fluid-constants", "refGrav", refGrav);
  r_ynolds = refRho * refVel * refLen / refMu;
  weber = refRho * refVel * refVel * refLen / refSurfTen;
  froude = refVel / sqrt(refGrav * refLen);

  // Temperature stuff
  coeff_thermal_expan = 3.43e-3;
  config->getDouble("fluid-constants", "coeff_thermal_expan", coeff_thermal_expan);
  DG_FP thermal_conductivity = 0.026; // W / mK
  config->getDouble("fluid-constants", "thermal_conductivity", thermal_conductivity);
  DG_FP heat_capacity = 0.000280 * 3.6e6; // J/kgK
  config->getDouble("fluid-constants", "heat_capacity", heat_capacity);

  // TODO: Check if this should include refRho
  // peclet = (heat_capacity * refVel * refLen) / thermal_conductivity;
  peclet = (heat_capacity * refVel * refLen * refRho) / thermal_conductivity;

  int re = -1;
  PetscOptionsGetInt(NULL, NULL, "-re", &re, &found);
  if(re > 0) {
    r_ynolds = (DG_FP)re;
  }
  op_printf("\n\nRe: %g\n", r_ynolds);
  op_printf("Weber: %g\n", weber);
  op_printf("Froude: %g\n", froude);
  op_printf("Peclet: %g\n\n", peclet);

  // For 2D compressible Euler
  gamma_e = 1.4;

  DGMesh2D *mesh = new DGMesh2D(filename);

  std::string type_of_solver = "single-phase";
  config->getStr("simulation-constants", "type_of_solver", type_of_solver);
  SimulationDriver* driver;
  vector<Measurement2D*> measurements;
  if(type_of_solver == "single-phase") {
    INSSolverBase2D *ins2d;
    if(resuming)
      ins2d = new INSSolver2D(mesh, r_ynolds, checkpointFile, resumeIter);
    else
      ins2d = new INSSolver2D(mesh, r_ynolds);
    driver = ins2d;
  } else if(type_of_solver == "multi-phase") {
    INSSolverBase2D *ins2d;
    if(resuming)
      ins2d = new MPINSSolver2D(mesh, r_ynolds, checkpointFile, resumeIter);
    else
      ins2d = new MPINSSolver2D(mesh, r_ynolds);
    driver = ins2d;
  } else if(type_of_solver == "single-phase-with-temperature") {
    INSSolverBase2D *ins2d;
    if(resuming)
      ins2d = new INSTemperatureSolver2D(mesh, r_ynolds, checkpointFile, resumeIter);
    else
      ins2d = new INSTemperatureSolver2D(mesh, r_ynolds);
    driver = ins2d;
  } else if(type_of_solver == "level-set-only") {
    if(resuming)
      throw std::runtime_error("Resuming level-set-only simulation is not supported currently");
    else
      driver = new LSDriver2D(mesh);
  } else if(type_of_solver == "compressible-euler") {
    if(resuming)
      throw std::runtime_error("Resuming level-set-only simulation is not supported currently");
    else
      driver = new CompressibleEulerDriver2D(mesh);
  } else {
    dg_abort("Unknown \'type_of_solver\' specified in config file, valid options: \'single-phase\', \'multi-phase\', \'single-phase-with-temperature\'");
  }

  add_measurements(driver, measurements, refMu);

  // Toolkit constants
  op_decl_const(DG_ORDER * 5, "int", DG_CONSTANTS);
  op_decl_const(DG_ORDER * 3 * DG_NPF, "int", FMASK);

  // Application constants
  op_decl_const(1, DG_FP_STR, &r_ynolds);
  op_decl_const(1, DG_FP_STR, &mu0);
  op_decl_const(1, DG_FP_STR, &mu1);
  op_decl_const(1, DG_FP_STR, &rho0);
  op_decl_const(1, DG_FP_STR, &rho1);
  op_decl_const(1, DG_FP_STR, &gamma_e);
  op_decl_const(1, DG_FP_STR, &weber);
  op_decl_const(1, DG_FP_STR, &froude);
  op_decl_const(1, DG_FP_STR, &coeff_thermal_expan);
  op_decl_const(1, DG_FP_STR, &peclet);

  timer->startTimer("OP2 Partitioning");
  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->face2cells, NULL);
  int renumber_elements = 0;
  config->getInt("simulation-constants", "renumber_elements", renumber_elements);
  if(renumber_elements) {
    op_renumber(mesh->face2cells);
  }
  timer->endTimer("OP2 Partitioning");

  mesh->init();
  driver->init();

  int save_ic = 1;
  config->getInt("io", "save_ic", save_ic);
  if(save_ic) {
    string out_file_ic = outputDir + "ic.h5";
    driver->dump_visualisation_data(out_file_ic);
  }

  timer->startTimer("Main loop");
  int i = 0;
  if(resuming) i = resumeIter;
  for(; i < iter; i++) {
    driver->step();

    if(save > 0 && (i + 1) % save == 0) {
      string out_file_tmp = outputDir + "iter-" + to_string(i + 1) + ".h5";
      driver->dump_visualisation_data(out_file_tmp);
    }

    if(checkpoint > 0 && (i + 1) % checkpoint == 0) {
      string out_file_tmp = outputDir + "checkpoint-" + to_string(i + 1) + ".h5";
      driver->dump_checkpoint_data(out_file_tmp);
      int tmp_iter = i + 1;
      op_write_const_hdf5("iter", 1, "int", (char *)&tmp_iter, out_file_tmp.c_str());
    }

    for(auto &measurement : measurements) {
      measurement->measure();
    }
  }
  timer->endTimer("Main loop");

  int save_end = 1;
  config->getInt("io", "save_end", save_end);
  if(save_end) {
    string out_file_end = outputDir + "end.h5";
    driver->dump_visualisation_data(out_file_end);
  }

  for(auto &measurement : measurements) {
    measurement->output(outputDir);
  }

  timer->exportTimings(outputDir + "timings", iter, driver->get_time());

  for(auto &measurement : measurements) {
    delete measurement;
  }

  // Print closing summary
  op_printf("\n\n Summary of simulation:\n");
  op_printf("%d iterations\n", iter);
  op_printf("%g time (non-dimensionalised)\n", driver->get_time());
  op_printf("%g time (s)\n", driver->get_time() * refLen / refVel);
  op_printf("Reference density: %g kg m^-3\n", refRho);
  op_printf("Reference velocity: %g m s^-1\n", refVel);
  op_printf("Reference length: %g m\n", refLen);
  op_printf("Reference viscosity: %g m^2 s^-1\n", refMu);
  op_printf("Reference gravity: %g m s^-2\n", refMu);
  op_printf("Reference surface tension: %g N m^-1\n", refSurfTen);
  op_printf("Density ratio of %g : %g\n", rho0, rho1);
  op_printf("Viscosity ratio of %g : %g\n", mu0, mu1);
  op_printf("Re (fluid 0): %g\n", r_ynolds * (rho0 / mu0));
  op_printf("Re (fluid 1): %g\n", r_ynolds * (rho1 / mu1));
  op_printf("Weber: %g\n", weber);
  op_printf("Froude: %g\n", froude);
  op_printf("Peclet: %g\n\n", peclet);

  delete driver;

  ierr = PetscFinalize();

  delete timer;
  delete config;

  // Clean up OP2
  op_exit();
  return ierr;
}
