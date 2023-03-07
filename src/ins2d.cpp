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
#include "solvers/2d/mp_ins_solver.h"
#include "solvers/2d/ins_solver.h"

Timing *timer;

using namespace std;

// Global constants
DG_FP r_ynolds, mu0, mu1, rho0, rho1, gamma_e;

int main(int argc, char **argv) {
  op_init(argc, argv, 1);

  timer = new Timing();

  char help[] = "Run for i iterations with \"-iter i\"\nSave solution every x iterations with \"-save x\"\n";
  int ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if(ierr) {
    op_printf("Error initialising PETSc\n");
    return ierr;
  }

  // Get input from args
  int iter = 1;
  PetscBool found;
  PetscOptionsGetInt(NULL, NULL, "-iter", &iter, &found);

  int save = -1;
  PetscOptionsGetInt(NULL, NULL, "-save", &save, &found);

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

  mu0  = 1.0;
  mu1  = 1.0;
  rho0 = 1.0;
  rho1 = 1.0;

  const DG_FP refRho = 1.0;
  const DG_FP refVel = 1.0;
  const DG_FP refLen = 0.005;
  const DG_FP refMu  = 1.0e-5;
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
  INSSolver2D *mpins2d;
  if(resumeIter == 0)
    mpins2d = new INSSolver2D(mesh);
  else
    mpins2d = new INSSolver2D(mesh, checkpointFile, resumeIter);

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

    if((i + 1) % 25 == 0) {
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

  // Clean up OP2
  op_exit();
  return ierr;
}
