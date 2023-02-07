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
#include "mp_ins_solver.h"

Timing *timer;

using namespace std;

// Global constants
double r_ynolds, mu0, mu1, rho0, rho1;

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

  const double refRho = 1.0;
  const double refVel = 1.0;
  const double refLen = 0.005;
  const double refMu  = 1.0e-5;
  r_ynolds = refRho * refVel * refLen / refMu;

  DGMesh2D *mesh = new DGMesh2D(filename);
  MPINSSolver2D *mpins2d = new MPINSSolver2D(mesh);

  // Toolkit constants
  op_decl_const(DG_ORDER * 5, "int", DG_CONSTANTS);
  op_decl_const(DG_ORDER * 3 * DG_NPF, "int", FMASK);
  op_decl_const(DG_ORDER * DG_CUB_NP, "double", cubW_g);
  op_decl_const(DG_ORDER * DG_GF_NP, "double", gaussW_g);

  // Application constants
  op_decl_const(1, "double", &r_ynolds);
  op_decl_const(1, "double", &mu0);
  op_decl_const(1, "double", &mu1);
  op_decl_const(1, "double", &rho0);
  op_decl_const(1, "double", &rho1);

  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->face2cells, NULL);

  mesh->init();
  mpins2d->init(r_ynolds, refVel);

  string out_file_ic = outputDir + "ls_test_ic.h5";
  mpins2d->dump_data(out_file_ic);

  for(int i = 0; i < iter; i++) {
    mpins2d->step();

    op_printf("Iter %d\n", i);
  }

  string out_file_end = outputDir + "ls_test_end.h5";
  mpins2d->dump_data(out_file_end);

  timer->exportTimings(outputDir + "timings.txt", 0, 0.0);

  delete mpins2d;

  ierr = PetscFinalize();

  delete timer;

  // Clean up OP2
  op_exit();
  return ierr;
}
