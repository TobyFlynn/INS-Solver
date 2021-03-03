// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

#include "constants/all_constants.h"
#include "load_mesh.h"
#include "ins_data.h"
#include "blas_calls.h"
#include "operators.h"
#include "save_solution.h"
#include "poisson.h"

// Kernels
#include "kernels/init_grid.h"
#include "kernels/poisson_test_init.h"
#include "kernels/poisson_test_bc.h"
#include "kernels/poisson_test_set_rhs.h"
#include "kernels/poisson_test_error.h"

using namespace std;

int main(int argc, char **argv) {
  char help[] = "TODO";
  int ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if(ierr) {
    cout << "Error initialising PETSc" << endl;
    return ierr;
  }

  // Object that holds all sets, maps and dats
  // (along with memory associated with them)
  INSData *data = new INSData();

  auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
    if(y1 == y2 && y1 > 0.5) {
      // Neumann BC y = 1
      // cout << "0" << endl;
      return 2;
    } else if(y1 == y2 && y1 < 0.5) {
      // Neumann BC y = 0
      // cout << "1" << endl;
      return 3;
    } else if(x1 < 0.5){
      // Dirichlet BC x = 0
      // cout << "2" << endl;
      return 0;
    } else {
      // Dirichlet BC x = 1
      // cout << "3" << endl;
      return 1;
    }
  };

  load_mesh("./grid.cgns", data, bcNum);

  // Initialise OP2
  op_init(argc, argv, 2);

  // Initialise all sets, maps and dats
  data->initOP2();

  double *rhs_data = (double *)malloc(15 * data->numCells * sizeof(double));
  double *sol_data = (double *)malloc(15 * data->numCells * sizeof(double));
  double *ex_data  = (double *)malloc(15 * data->numCells * sizeof(double));
  double *err_data = (double *)malloc(15 * data->numCells * sizeof(double));

  op_dat rhs = op_decl_dat(data->cells, 15, "double", rhs_data, "rhs");
  op_dat sol = op_decl_dat(data->cells, 15, "double", sol_data, "sol");
  op_dat ex  = op_decl_dat(data->cells, 15, "double", ex_data, "ex");
  op_dat err = op_decl_dat(data->cells, 15, "double", err_data, "err");

  // Calculate geometric factors
  init_grid_blas(data);

  op_par_loop(init_grid, "init_grid", data->cells,
              op_arg_dat(data->node_coords, -3, data->cell2nodes, 2, "double", OP_READ),
              op_arg_dat(data->nodeX, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(data->nodeY, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(data->xr, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->yr, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->xs, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ys, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->rx, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->ry, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->sx, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->sy, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->nx, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->ny, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->J,  -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->sJ, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->fscale, -1, OP_ID, 15, "double", OP_WRITE));

  op_par_loop(poisson_test_init, "poisson_test_init", data->cells,
              op_arg_dat(data->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(ex,  -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->dirichletBC, -1, OP_ID, 15, "double", OP_WRITE));

  op_par_loop(poisson_test_bc, "poisson_test_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->x, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->y, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->dirichletBC, 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(poisson_test_set_rhs, "poisson_test_set_rhs", data->cells,
              op_arg_dat(data->J,  -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(ex,  -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_RW));

  poisson_test_rhs_blas(data, rhs);

  Poisson *poisson = new Poisson(data);

  // op_fetch_data_hdf5_file(data->dirichletBC, "t.h5");

  int dBCs[] = {0, 1};
  int nBCs[] = {2, 3};
  poisson->setDirichletBCs(dBCs, data->dirichletBC);
  poisson->setNeumannBCs(nBCs);

  poisson->solve(rhs, sol, false);

  double l2 = 0.0;
  op_par_loop(poisson_test_error, "poisson_test_error", data->cells,
              op_arg_dat(data->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(sol, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(err, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_gbl(&l2, 1, "double", OP_INC));

  cout << "L^2 error: " << l2 << endl;

  // Save solution to CGNS file (get it twice just to avoid rewriting save_solution for this test)
  double *sol_ptr = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *rhs_ptr = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *err_ptr = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *x_ptr = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  double *y_ptr = (double *)malloc(15 * op_get_size(data->cells) * sizeof(double));
  op_fetch_data(sol, sol_ptr);
  op_fetch_data(rhs, rhs_ptr);
  op_fetch_data(err, err_ptr);
  op_fetch_data(data->x, x_ptr);
  op_fetch_data(data->y, y_ptr);
  save_solution_all("./sol.cgns", op_get_size(data->cells), sol_ptr, err_ptr, x_ptr, y_ptr);

  save_solution("./grid.cgns", op_get_size(data->nodes), op_get_size(data->cells),
                sol_ptr, err_ptr, rhs_ptr, data->cgnsCells);

  // save_solution_cell("./grid.cgns", op_get_size(data->nodes), op_get_size(data->cells),
  //               sol_ptr, err_ptr, data->cgnsCells);

  free(sol_ptr);
  free(rhs_ptr);
  free(err_ptr);
  free(x_ptr);
  free(y_ptr);

  // Clean up OP2
  op_exit();

  delete poisson;
  delete data;

  free(rhs_data);
  free(sol_data);
  free(ex_data);
  free(err_data);

  ierr = PetscFinalize();
  return ierr;
}
