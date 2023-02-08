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

  DGMesh3D *mesh = new DGMesh3D(filename);
  AdvectionSolver3D *advec3d = new AdvectionSolver3D(mesh);

  double *data_t0 = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  op_dat u   = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advec_u");
  op_dat v   = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advec_v");
  op_dat w   = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advec_w");
  op_dat val = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advec_val");
  free(data_t0);

  // Toolkit constants
  op_decl_const(DG_ORDER * 2, "int", DG_CONSTANTS);
  op_decl_const(DG_ORDER * 4 * DG_NPF, "int", FMASK);

  // Application constants
  op_decl_const(1, "double", &r_ynolds);
  op_decl_const(1, "double", &mu0);
  op_decl_const(1, "double", &mu1);
  op_decl_const(1, "double", &rho0);
  op_decl_const(1, "double", &rho1);

  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->face2cells, NULL);

  mesh->init();

  op_par_loop(advec_test_3d, "advec_test_3d", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->z, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u,   -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(v,   -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(w,   -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(val, -1, OP_ID, DG_NP, "double", OP_WRITE));

  string out_file_ic = outputDir + "test_ic.h5";
  op_fetch_data_hdf5_file(mesh->x, out_file_ic.c_str());
  op_fetch_data_hdf5_file(mesh->y, out_file_ic.c_str());
  op_fetch_data_hdf5_file(mesh->z, out_file_ic.c_str());
  op_fetch_data_hdf5_file(u, out_file_ic.c_str());
  op_fetch_data_hdf5_file(v, out_file_ic.c_str());
  op_fetch_data_hdf5_file(w, out_file_ic.c_str());
  op_fetch_data_hdf5_file(val, out_file_ic.c_str());

  for(int i = 0; i < iter; i++) {
    advec3d->step(val, u, v, w);

    op_printf("Iter %d\n", i);
  }

  string out_file_end = outputDir + "test_end.h5";
  op_fetch_data_hdf5_file(mesh->x, out_file_end.c_str());
  op_fetch_data_hdf5_file(mesh->y, out_file_end.c_str());
  op_fetch_data_hdf5_file(mesh->z, out_file_end.c_str());
  op_fetch_data_hdf5_file(u, out_file_end.c_str());
  op_fetch_data_hdf5_file(v, out_file_end.c_str());
  op_fetch_data_hdf5_file(w, out_file_end.c_str());
  op_fetch_data_hdf5_file(val, out_file_end.c_str());


  timer->exportTimings(outputDir + "timings.txt", 0, 0.0);

  delete advec3d;

  ierr = PetscFinalize();

  delete timer;

  // Clean up OP2
  op_exit();
  return ierr;
}
