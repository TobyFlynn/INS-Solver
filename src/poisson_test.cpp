#include "op_seq.h"

#include <iostream>

#include "dg_constants.h"
#include "dg_mesh.h"
#include "ins_data.h"
#include "save_solution.h"
#include "poisson.h"
#include "timing.h"
#include "dg_blas_calls.h"
#include "load_mesh.h"
#include "dg_global_constants.h"

using namespace std;

Timing *timer;
DGConstants *constants;

int main(int argc, char **argv) {
  // Initialise OP2
  op_init(argc, argv, 2);

  timer = new Timing();
  timer->startWallTime();
  timer->startSetup();
  constants = new DGConstants();

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

  double *coords_data;
  int *cells_data, *edge2node_data, *edge2cell_data, *bedge2node_data;
  int *bedge2cell_data, *bedge_type_data, *edgeNum_data, *bedgeNum_data;
  int numNodes_g, numCells_g, numEdges_g, numBoundaryEdges_g, numNodes;
  int numCells, numEdges, numBoundaryEdges;

  int pressure_dirichlet[3];
  int pressure_neumann[3];
  int viscosity_dirichlet[3];
  int viscosity_neumann[3];

  load_mesh(filename, &coords_data, &cells_data, &edge2node_data,
            &edge2cell_data, &bedge2node_data, &bedge2cell_data,
            &bedge_type_data, &edgeNum_data, &bedgeNum_data, &numNodes_g,
            &numCells_g, &numEdges_g, &numBoundaryEdges_g, &numNodes, &numCells,
            &numEdges, &numBoundaryEdges, pressure_dirichlet, pressure_neumann,
            viscosity_dirichlet, viscosity_neumann);

  DGMesh *mesh = new DGMesh(coords_data, cells_data, edge2node_data, edge2cell_data,
                    bedge2node_data, bedge2cell_data, bedge_type_data,
                    edgeNum_data, bedgeNum_data, numNodes_g, numCells_g,
                    numEdges_g, numBoundaryEdges_g, numNodes, numCells,
                    numEdges, numBoundaryEdges);
  INSData *data = new INSData(mesh);

  Poisson *poisson;
  int dirichlet[] = {0, 1, -1};
  int neumann[] = {2, 3, -1};

  if(pmethod == 0) {
    Poisson_M *pressureM = new Poisson_M(mesh, data);
    pressureM->setDirichletBCs(dirichlet);
    pressureM->setNeumannBCs(neumann);
    poisson = pressureM;
  } else if(pmethod == 1) {
    Poisson_MF *poissonMF = new Poisson_MF(mesh, data);
    poissonMF->setDirichletBCs(dirichlet);
    poissonMF->setNeumannBCs(neumann);
    poisson = poissonMF;
  } else {
    Poisson_MF2 *poissonMF2 = new Poisson_MF2(mesh, data);
    poissonMF2->setDirichletBCs(dirichlet);
    poissonMF2->setNeumannBCs(neumann);
    poisson = poissonMF2;
  }

  double *tmp_data = (double *)calloc(15 * mesh->numCells, sizeof(double));
  double *rhs_data = (double *)calloc(15 * mesh->numCells, sizeof(double));
  double *bc_data  = (double *)calloc(21 * mesh->numCells, sizeof(double));

  op_dat tmp = op_decl_dat(mesh->cells, 15, "double", tmp_data, "tmp");
  op_dat rhs = op_decl_dat(mesh->cells, 15, "double", rhs_data, "rhs");
  op_dat bc  = op_decl_dat(mesh->cells, 21, "double", bc_data, "bc");

  op_partition("PARMETIS", "KWAY", mesh->cells, mesh->edge2cells, NULL);

  mesh->init();
  data->init();
  poisson->init();

  op_par_loop(poisson_test_init, "poisson_test_init", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(tmp, -1, OP_ID, 15, "double", OP_WRITE));

  op_par_loop(poisson_test_bc, "poisson_test_bc", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(bc, 0, mesh->bedge2cells, 21, "double", OP_INC));

  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(DGConstants::MASS), 15, tmp, 0.0, rhs);

  poisson->setBCValues(bc);
  poisson->solve(rhs, data->p);

  op_par_loop(poisson_test_error, "poisson_test_error", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->p, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->Q[0][0], -1, OP_ID, 15, "double", OP_WRITE));

  save_solution("end.cgns", mesh, data, 0, NULL);

  free(rhs_data);
  free(bc_data);

  delete poisson;
  delete data;
  delete mesh;
  delete constants;
  delete timer;

  ierr = PetscFinalize();
  // Clean up OP2
  op_exit();
  return ierr;
}
