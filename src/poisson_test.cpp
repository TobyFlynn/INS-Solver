#include "op_seq.h"

#include <iostream>

#include "constants/all_constants.h"
#include "constants.h"
#include "ins_data.h"
#include "save_solution.h"
#include "poisson.h"
#include "timing.h"
#include "blas_calls.h"

#include "kernels/poisson_test_init.h"
#include "kernels/poisson_test_bc.h"
#include "kernels/poisson_test_error.h"
#include "kernels/poisson_test_tmp.h"

using namespace std;

Timing *timer;
Constants *constants;

int main(int argc, char **argv) {
  // Initialise OP2
  op_init(argc, argv, 2);

  timer = new Timing();
  timer->startWallTime();
  timer->startSetup();
  constants = new Constants();

  string filename = "./grid.cgns";
  char help[] = "TODO\n";
  int ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if(ierr) {
    cout << "Error initialising PETSc" << endl;
    return ierr;
  }

  PetscBool found;
  int pmethod = 0;
  PetscOptionsGetInt(NULL, NULL, "-pmethod", &pmethod, &found);

  INSData *data = new INSData(filename);
  CubatureData *cubData = new CubatureData(data);
  GaussData *gaussData = new GaussData(data);

  Poisson *poisson;
  int dirichlet[] = {0, 1, -1};
  int neumann[] = {2, 3, -1};

  if(pmethod == 0) {
    Poisson_M *pressureM = new Poisson_M(data, cubData, gaussData);
    pressureM->setDirichletBCs(dirichlet);
    pressureM->setNeumannBCs(neumann);
    poisson = pressureM;
  } else {
    Poisson_MF2 *poissonMF2 = new Poisson_MF2(data, cubData, gaussData);
    poissonMF2->setDirichletBCs(dirichlet);
    poissonMF2->setNeumannBCs(neumann);
    poisson = poissonMF2;
  }

  double *tmp_data = (double *)calloc(15 * data->numCells, sizeof(double));
  double *rhs_data = (double *)calloc(15 * data->numCells, sizeof(double));
  double *bc_data  = (double *)calloc(21 * data->numCells, sizeof(double));

  op_dat tmp = op_decl_dat(data->cells, 15, "double", tmp_data, "tmp");
  op_dat rhs = op_decl_dat(data->cells, 15, "double", rhs_data, "rhs");
  op_dat bc  = op_decl_dat(data->cells, 21, "double", bc_data, "bc");

  op_partition("PARMETIS", "KWAY", data->cells, data->edge2cells, NULL);
  // op_partition("PARMETIS", "KWAY", data->cells, NULL, NULL);

  data->init();
  cubData->init();
  gaussData->init();
  poisson->init();

  op_par_loop(poisson_test_init, "poisson_test_init", data->cells,
              op_arg_dat(data->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(tmp, -1, OP_ID, 15, "double", OP_WRITE));

  op_par_loop(poisson_test_bc, "poisson_test_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(gaussData->x, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gaussData->y, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(bc, 0, data->bedge2cells, 21, "double", OP_INC));

  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::MASS), 15, tmp, 0.0, rhs);

  poisson->setBCValues(bc);
  poisson->solve(rhs, data->p);

  op_par_loop(poisson_test_error, "poisson_test_error", data->cells,
              op_arg_dat(data->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->p, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[0][0], -1, OP_ID, 15, "double", OP_WRITE));

  save_solution("end.cgns", data, 0, NULL);

  free(rhs_data);
  free(bc_data);

  delete poisson;
  delete gaussData;
  delete cubData;
  delete data;
  delete constants;
  delete timer;

  ierr = PetscFinalize();
  // Clean up OP2
  op_exit();
  return ierr;
}
