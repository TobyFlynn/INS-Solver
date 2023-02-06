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
#include "advection_solver.h"
#include "ls_solver.h"

Timing *timer;

using namespace std;

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

  DGMesh2D *mesh = new DGMesh2D(filename);
  // AdvectionSolver2D *advec2d = new AdvectionSolver2D(mesh);
  LevelSetSolver2D *ls2d = new LevelSetSolver2D(mesh);

  double *data_t0 = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  op_dat u = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "ls_test_u");
  op_dat v = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "ls_test_v");
  free(data_t0);

  op_decl_const(DG_ORDER * 5, "int", DG_CONSTANTS);
  op_decl_const(DG_ORDER * 3 * DG_NPF, "int", FMASK);
  op_decl_const(DG_ORDER * DG_CUB_NP, "double", cubW_g);
  op_decl_const(DG_ORDER * DG_GF_NP, "double", gaussW_g);

  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->face2cells, NULL);

  mesh->init();
  ls2d->init();

  op_par_loop(ls_test_ic, "ls_test_ic", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u,   -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(v,   -1, OP_ID, DG_NP, "double", OP_WRITE));

  string out_file_ic = outputDir + "ls_test_ic.h5";
  op_fetch_data_hdf5_file(mesh->x, out_file_ic.c_str());
  op_fetch_data_hdf5_file(mesh->y, out_file_ic.c_str());
  op_fetch_data_hdf5_file(ls2d->s, out_file_ic.c_str());
  op_fetch_data_hdf5_file(u, out_file_ic.c_str());
  op_fetch_data_hdf5_file(v, out_file_ic.c_str());

  for(int i = 0; i < iter; i++) {
    ls2d->setVelField(u, v);
    ls2d->step(1e-3);

    op_printf("Iter %d\n", i);
  }

  string out_file_end = outputDir + "ls_test_end.h5";
  op_fetch_data_hdf5_file(mesh->x, out_file_end.c_str());
  op_fetch_data_hdf5_file(mesh->y, out_file_end.c_str());
  op_fetch_data_hdf5_file(ls2d->s, out_file_end.c_str());
  op_fetch_data_hdf5_file(u, out_file_end.c_str());
  op_fetch_data_hdf5_file(v, out_file_end.c_str());

  timer->exportTimings(outputDir + "timings.txt", 0, 0.0);

  ierr = PetscFinalize();

  delete timer;

  // Clean up OP2
  op_exit();
  return ierr;
}
