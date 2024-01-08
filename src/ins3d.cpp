#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>

#include "petscvec.h"
#include "petscksp.h"

#include "dg_global_constants/dg_global_constants_3d.h"
#include "dg_dat_pool.h"
#include "dg_abort.h"

#include "timing.h"
#include "config.h"
#include "op2_utils.h"

#include "solvers/3d/advection_solver.h"
#include "solvers/3d/ins_solver.h"
#include "solvers/3d/mp_ins_solver.h"
#include "measurements/3d/enstropy.h"

Timing *timer;
Config *config;
extern DGDatPool *dg_dat_pool;

using namespace std;

// Global constants
DG_FP r_ynolds, mu0, mu1, rho0, rho1;

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

  // Resuming from checkpoint
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

  // Get input from args
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
  DG_FP refLen = 0.001;
  DG_FP refMu  = 1.0e-5;
  config->getDouble("fluid-constants", "refRho", refRho);
  config->getDouble("fluid-constants", "refVel", refVel);
  config->getDouble("fluid-constants", "refLen", refLen);
  config->getDouble("fluid-constants", "refMu", refMu);
  r_ynolds = refRho * refVel * refLen / refMu;

  int re = -1;
  PetscOptionsGetInt(NULL, NULL, "-re", &re, &found);
  if(re > 0) {
    r_ynolds = (DG_FP)re;
  }

  op_printf("Reynolds number: %g\n", r_ynolds);

  DGMesh3D *mesh = new DGMesh3D(filename);

  std::string type_of_solver = "single-phase";
  config->getStr("simulation-constants", "type_of_solver", type_of_solver);
  INSSolverBase3D *ins3d;
  if(type_of_solver == "single-phase") {
    if(resuming)
      ins3d = new INSSolver3D(mesh, checkpointFile, resumeIter);
    else
      ins3d = new INSSolver3D(mesh);
  } else if(type_of_solver == "multi-phase") {
    if(resuming)
      ins3d = new MPINSSolver3D(mesh, checkpointFile, resumeIter);
    else
      ins3d = new MPINSSolver3D(mesh);
  } else {
    dg_abort("Unknown \'type_of_solver\' specified in config file, valid options: \'single-phase\', \'multi-phase\'");
  }

  // Toolkit constants
  op_decl_const(DG_ORDER * 2, "int", DG_CONSTANTS);
  op_decl_const(DG_ORDER * 4 * DG_NPF, "int", FMASK);

  // Application constants
  op_decl_const(1, DG_FP_STR, &r_ynolds);
  op_decl_const(1, DG_FP_STR, &mu0);
  op_decl_const(1, DG_FP_STR, &mu1);
  op_decl_const(1, DG_FP_STR, &rho0);
  op_decl_const(1, DG_FP_STR, &rho1);

  timer->startTimer("OP2 Partitioning");
  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->face2cells, NULL);
  int renumber_elements = 0;
  config->getInt("simulation-constants", "renumber_elements", renumber_elements);
  if(renumber_elements) {
    op_renumber(mesh->face2cells);
  }
  timer->endTimer("OP2 Partitioning");

  mesh->init();
  ins3d->init(r_ynolds, refVel);

  // Set up measurements
  vector<Measurement3D*> measurements;
  // Get list of measurements to take
  // Options are: enstropy
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
      if(measurement == "enstropy") {
        Enstropy3D *enstropy = new Enstropy3D(ins3d, refMu, refRho, refVel, 248.0502134423985614038105205);
        measurements.push_back(enstropy);
      } else {
        dg_abort("Unrecognised measurement: " + measurement);
      }
    }
  }

  int save_ic = 1;
  config->getInt("io", "save_ic", save_ic);
  if(save_ic) {
    string out_file_ic = outputDir + "ic.h5";
    ins3d->dump_visualisation_data(out_file_ic);
  }

  timer->startTimer("Main loop");
  int i = 0;
  if(resuming) i = resumeIter;
  for(; i < iter; i++) {
    ins3d->step();

    if(save > 0 && (i + 1) % save == 0) {
      string out_file_tmp = outputDir + "iter-" + to_string(i + 1) + ".h5";
      ins3d->dump_visualisation_data(out_file_tmp);
    }

    if(checkpoint > 0 && (i + 1) % checkpoint == 0) {
      string out_file_tmp = outputDir + "checkpoint-" + to_string(i + 1) + ".h5";
      ins3d->dump_checkpoint_data(out_file_tmp);
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
    ins3d->dump_visualisation_data(out_file_end);
  }

  for(auto &measurement : measurements) {
    measurement->output(outputDir);
  }

  timer->exportTimings(outputDir + "timings", iter, ins3d->get_time());

  for(auto &measurement : measurements) {
    delete measurement;
  }

  dg_dat_pool->report();

  // Print closing summary
  op_printf("\n\n Summary of simulation:\n");
  op_printf("%d iterations\n", iter);
  op_printf("%g time (non-dimensionalised)\n", ins3d->get_time());
  op_printf("%g time (s)\n", ins3d->get_time() * refLen / refVel);
  op_printf("Reference density: %g kg m^-3\n", refRho);
  op_printf("Reference velocity: %g m s^-1\n", refVel);
  op_printf("Reference length: %g m\n", refLen);
  op_printf("Reference viscosity: %g m^2 s^-1\n", refMu);
  op_printf("Density ratio of %g : %g\n", rho0, rho1);
  op_printf("Viscosity ratio of %g : %g\n", mu0, mu1);
  op_printf("Re: %g\n", r_ynolds);

  delete ins3d;

  ierr = PetscFinalize();

  delete config;
  delete timer;

  // Clean up OP2
  op_exit();
  return ierr;
}
