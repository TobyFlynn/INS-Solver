#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>

#include "petscvec.h"
#include "petscksp.h"

#include "dg_global_constants/dg_global_constants_3d.h"

#include "timing.h"
#include "solvers/3d/advection_solver.h"
#include "solvers/3d/ins_solver.h"

Timing *timer;

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

  mu0  = 1.0;
  mu1  = 1.0;
  rho0 = 1.0;
  rho1 = 1.0;

  const DG_FP refRho = 1.0;
  const DG_FP refVel = 1.0;
  const DG_FP refLen = 0.005;
  const DG_FP refMu  = 1.0e-5;
  r_ynolds = refRho * refVel * refLen / refMu;

  op_printf("Reynolds number: %g\n", r_ynolds);

  DGMesh3D *mesh = new DGMesh3D(filename);
  INSSolver3D *ins3d = new INSSolver3D(mesh);

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
  timer->endTimer("OP2 Partitioning");

  mesh->init();
  ins3d->init(r_ynolds, refVel);

  string out_file_ic = outputDir + "test_ic.h5";
  ins3d->dump_data(out_file_ic);

  timer->startTimer("Main loop");
  for(int i = 0; i < iter; i++) {
    ins3d->step();
  }
  timer->endTimer("Main loop");

  string out_file_end = outputDir + "test_end.h5";
  ins3d->dump_data(out_file_end);


  timer->exportTimings(outputDir + "timings.txt", iter, ins3d->get_time());

  delete ins3d;

  ierr = PetscFinalize();

  delete timer;

  // Clean up OP2
  op_exit();
  return ierr;
}