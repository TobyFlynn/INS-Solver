#include "op_seq.h"

#include <iostream>

#include "constants/all_constants.h"
#include "constants.h"
#include "load_mesh.h"
#include "ins_data.h"
#include "save_solution.h"
#include "poisson.h"
#include "blas_calls.h"
#include "timing.h"

#include "kernels/poisson_test_init.h"
#include "kernels/poisson_test_bc.h"
#include "kernels/poisson_test_error.h"

using namespace std;

Timing *timer;
Constants *constants;

int main(int argc, char **argv) {
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

  INSData *data = new INSData();

  // // Lamda used to identify the type of boundary edges
  // auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
  //   if(y1 == y2 && y1 > 0.5) {
  //     // Neumann BC y = 1
  //     return 2;
  //   } else if(y1 == y2 && y1 < 0.5) {
  //     // Neumann BC y = 0
  //     return 3;
  //   } else if(x1 < 0.5){
  //     // Dirichlet BC x = 0
  //     return 0;
  //   } else {
  //     // Dirichlet BC x = 1
  //     return 1;
  //   }
  // };

  // Lamda used to identify the type of boundary edges
  auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
    if(y1 == y2 && y1 > 0.5) {
      // Neumann BC y = 1
      return 1;
    } else if(y1 == y2 && y1 < 0.5) {
      // Neumann BC y = 0
      return 1;
    } else if(x1 < 0.5){
      // Dirichlet BC x = 0
      return 1;
    } else {
      // Dirichlet BC x = 1
      return 0;
    }
  };

  load_mesh(filename.c_str(), data, bcNum);

  // Initialise all sets, maps and dats
  data->initOP2();

  CubatureData *cubData = new CubatureData(data);
  GaussData *gaussData = new GaussData(data);

  double *rhs_data = (double *)calloc(15 * data->numCells, sizeof(double));
  double *bc_data  = (double *)calloc(21 * data->numCells, sizeof(double));

  op_dat rhs = op_decl_dat(data->cells, 15, "double", rhs_data, "rhs");
  op_dat bc  = op_decl_dat(data->cells, 21, "double", bc_data, "bc");

  op_par_loop(poisson_test_init, "poisson_test_init", data->cells,
              op_arg_dat(data->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_WRITE));

  op_par_loop(poisson_test_bc, "poisson_test_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(gaussData->x, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(gaussData->y, 0, data->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(bc, 0, data->bedge2cells, 21, "double", OP_INC));

  poisson_test_rhs_blas(data, rhs);

  Poisson_MF *poisson = new Poisson_MF(data, cubData, gaussData);
  int dirichlet[] = {0, 1, -1};
  int neumann[] = {2, 3, -1};
  poisson->setDirichletBCs(dirichlet);
  poisson->setNeumannBCs(neumann);
  poisson->createBCMatrix();
  poisson->setBCValues(bc);
  // Poisson_M *poisson = new Poisson_M(data, cubData, gaussData);
  // int dirichlet[] = {0, 1, -1};
  // int neumann[] = {2, 3, -1};
  // poisson->setDirichletBCs(dirichlet);
  // poisson->setNeumannBCs(neumann);
  // poisson->createMatrix();
  // poisson->createBCMatrix();
  // poisson->setBCValues(bc);

  poisson->solve(rhs, data->p);

  op_par_loop(poisson_test_error, "poisson_test_error", data->cells,
              op_arg_dat(data->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->p, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[0][0], -1, OP_ID, 15, "double", OP_WRITE));

  save_solution("end.cgns", data, 0);

  free(rhs_data);
  free(bc_data);

  delete poisson;
  delete gaussData;
  delete cubData;
  delete data;
  delete constants;
  delete timer;

  ierr = PetscFinalize();
  return ierr;
}