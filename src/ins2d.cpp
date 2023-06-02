#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>

#include "petscvec.h"
#include "petscksp.h"

#include "dg_global_constants/dg_global_constants_2d.h"

#include "timing.h"
#include "config.h"
#include "solvers/2d/mp_ins_solver_over_int.h"
#include "solvers/2d/ins_solver_over_int.h"

Timing *timer;
Config *config;

using namespace std;

// Global constants
DG_FP r_ynolds, mu0, mu1, rho0, rho1, gamma_e;

int main(int argc, char **argv) {
  #ifdef DG_OP2_SOA
  throw std::runtime_error("2D not implemented for SoA");
  #endif

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

  int resumeIter = 0;
  PetscOptionsGetInt(NULL, NULL, "-r_iter", &resumeIter, &found);
  string checkpointFile;
  if(found) {
    // Resuming from checkpoint
    char checkFile[255];
    PetscOptionsGetString(NULL, NULL, "-checkpoint", checkFile, 255, &found);
    if(!found) {
      op_printf("Did not specify a checkpoint file after specifying iteration to resume from, use the -checkpoint flag\n");
      return -1;
    }
    checkpointFile = string(checkFile);
  }

  int iter = 1;
  config->getInt("simulation-constants", "iter", iter);

  int save = -1;
  config->getInt("simulation-constants", "save", save);

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
  op_printf("Re: %g\n", r_ynolds);

  // For 2D compressible Euler
  gamma_e = 1.4;

  DGMesh2D *mesh = new DGMesh2D(filename);
  INSSolverOverInt2D *mpins2d;
  if(resumeIter == 0)
    mpins2d = new INSSolverOverInt2D(mesh);
  else
    mpins2d = new INSSolverOverInt2D(mesh, checkpointFile, resumeIter);

  // Toolkit constants
  op_decl_const(DG_ORDER * 5, "int", DG_CONSTANTS);
  op_decl_const(DG_ORDER * 3 * DG_NPF, "int", FMASK);
  op_decl_const(DG_ORDER * DG_CUB_NP, DG_FP_STR, cubW_g);
  op_decl_const(DG_ORDER * DG_GF_NP, DG_FP_STR, gaussW_g);

  // Application constants
  op_decl_const(1, DG_FP_STR, &r_ynolds);
  op_decl_const(1, DG_FP_STR, &mu0);
  op_decl_const(1, DG_FP_STR, &mu1);
  op_decl_const(1, DG_FP_STR, &rho0);
  op_decl_const(1, DG_FP_STR, &rho1);
  op_decl_const(1, DG_FP_STR, &gamma_e);

  timer->startTimer("OP2 Partitioning");
  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->face2cells, NULL);
  timer->endTimer("OP2 Partitioning");

  mesh->init();
  mpins2d->init(r_ynolds, refVel);

  string out_file_ic = outputDir + "ic.h5";
  mpins2d->dump_data(out_file_ic);

  timer->startTimer("Main loop");
  for(int i = resumeIter; i < iter; i++) {
    mpins2d->step();

    op_printf("Iter %d\n", i);

    if(save > 0 && (i + 1) % save == 0) {
      string out_file_tmp = outputDir + "iter-" + to_string(i + 1) + ".h5";
      mpins2d->dump_data(out_file_tmp);
    }
  }
  timer->endTimer("Main loop");

  string out_file_end = outputDir + "end.h5";
  mpins2d->dump_data(out_file_end);

  timer->exportTimings(outputDir + "timings.txt", iter, mpins2d->get_time());

  delete mpins2d;

  ierr = PetscFinalize();

  delete timer;
  delete config;

  // Clean up OP2
  op_exit();
  return ierr;
}
