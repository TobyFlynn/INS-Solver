// Include OP2 stuff
#include "op_seq.h"
// Include CGNS stuff
// #include "cgnslib.h"
#include "pcgnslib.h"
#include "mpi.h"
// Include C++ stuff
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

#include "mpi_helper_func.h"
#include "dg_mesh.h"
#include "ins_data.h"
#include "ls.h"
#include "dg_operators.h"
#include "utils.h"
#include "shared_save_functions.h"

using namespace std;

void save_solution_iter(std::string filename, DGMesh *mesh, INSData *data, int ind, LS *ls, int iter) {
  // Calculate vorticity
  curl(mesh, data->Q[ind][0], data->Q[ind][1], data->vorticity);

  int file;
  if (cgp_open(filename.c_str(), CG_MODE_MODIFY, &file)) {
    cgp_error_exit();
  }
  int baseIndex = 1;
  int zoneIndex = 1;

  int flowIndex;
  string flowName = "FlowSolution" + to_string(iter);
  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, flowName.c_str(), CGNS_ENUMV(CellCenter), &flowIndex);

  cgsize_t numCells = mesh->cells->size * 9;
  cgsize_t minCell = 9 * get_global_start_index(mesh->cells) + 1;
  cgsize_t maxCell = minCell + numCells - 1;

  op_par_loop(save_values, "save_values", mesh->cells,
              op_arg_dat(data->Q[ind][0], -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

  int velXIndex;
  double *velX_data = getOP2Array(data->save_temp);
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", &velXIndex);
  cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, velXIndex, &minCell, &maxCell, velX_data);
  free(velX_data);

  op_par_loop(save_values, "save_values", mesh->cells,
              op_arg_dat(data->Q[ind][1], -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

  int velYIndex;
  double *velY_data = getOP2Array(data->save_temp);
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", &velYIndex);
  cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, velYIndex, &minCell, &maxCell, velY_data);
  free(velY_data);

  op_par_loop(save_values, "save_values", mesh->cells,
              op_arg_dat(data->p,         -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

  int pIndex;
  double *p_data = getOP2Array(data->save_temp);
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", &pIndex);
  cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, pIndex, &minCell, &maxCell, p_data);
  free(p_data);

  op_par_loop(save_values, "save_values", mesh->cells,
              op_arg_dat(data->vorticity, -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

  int vortIndex;
  double *vort_data = getOP2Array(data->save_temp);
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VorticityMagnitude", &vortIndex);
  cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, vortIndex, &minCell, &maxCell, vort_data);
  free(vort_data);

  int sIndex;
  if(ls) {
    op_par_loop(save_values, "save_values", mesh->cells,
                op_arg_dat(ls->s,           -1, OP_ID, 10, "double", OP_READ),
                op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

    double *s_data = getOP2Array(data->save_temp);
    cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Surface", &sIndex);
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, sIndex, &minCell, &maxCell, vort_data);
    free(s_data);
  }

  float exp[5];
  cg_goto(file, baseIndex, "end");
  cg_dataclass_write(CGNS_ENUMV(Dimensional));
  cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Degree));

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velXIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velYIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 1.0f; exp[1] = -1.0f; exp[2] = -2.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", pIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 0.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", vortIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  if(ls) {
    exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = 0.0f; exp[3] = 0.0f; exp[4] = 0.0f;
    cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", sIndex, "end");
    cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
  }

  cgp_close(file);
}

void save_solution_init(std::string filename, DGMesh *mesh, INSData *data, LS *ls) {
  // Calculate vorticity
  curl(mesh, data->Q[0][0], data->Q[0][1], data->vorticity);
  // Gather onto root process
  int rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  double *x_g;
  double *y_g;

  if(rank == 0) {
    x_g   = (double *)malloc(10 * mesh->numCells_g * sizeof(double));
    y_g   = (double *)malloc(10 * mesh->numCells_g * sizeof(double));
  }

  op_arg args[] = {
    op_arg_dat(mesh->x, -1, OP_ID, 10, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, 10, "double", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->cells, 2, args);

  gather_op2_double_array(x_g, (double *)mesh->x->data, mesh->cells->size, 10, comm_size, rank);
  gather_op2_double_array(y_g, (double *)mesh->y->data, mesh->cells->size, 10, comm_size, rank);

  op_mpi_set_dirtybit(2, args);

  vector<double> x_v;
  vector<double> y_v;
  vector<cgsize_t> cells(3 * mesh->numCells_g * 9);

  if(rank == 0) {
    if(DG_ORDER == 4) {
      get_cells_order_4(x_v, y_v, cells, x_g, y_g, mesh->numCells_g);
    } else if(DG_ORDER == 3) {
      get_cells_order_3(x_v, y_v, cells, x_g, y_g, mesh->numCells_g);
    } else if(DG_ORDER == 2) {
      get_cells_order_2(x_v, y_v, cells, x_g, y_g, mesh->numCells_g);
    } else {
      get_cells_order_1(x_v, y_v, cells, x_g, y_g, mesh->numCells_g);
    }
  }

  int file;
  if (cgp_open(filename.c_str(), CG_MODE_WRITE, &file)) {
    cgp_error_exit();
  }

  // Create base
  int baseIndex;
  int zoneIndex;
  int cellDim = 2;
  int physicalDim = 2;
  cg_base_write(file, "Base", cellDim, physicalDim, &baseIndex);

  int sizes_i[3];
  if(rank == 0) {
    sizes_i[0] = 10 * mesh->numCells_g;
    sizes_i[1] = mesh->numCells_g * 9;
    sizes_i[2] = 0;
  }
  MPI_Bcast(sizes_i, 3, MPI_INT, 0, MPI_COMM_WORLD);
  cgsize_t sizes[] = {sizes_i[0], sizes_i[1], sizes_i[2]};
  cg_zone_write(file, baseIndex, "Zone", sizes,
                CGNS_ENUMV(Unstructured), &zoneIndex);

  // Write grid coordinates
  int coordIndexX;
  int coordIndexY;
  cgp_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble), "CoordinateX", &coordIndexX);
  cgp_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble), "CoordinateY", &coordIndexY);

  if(rank == 0) {
    cgsize_t min = 1;
    cgsize_t max = x_v.size();
    cgp_coord_write_data(file, baseIndex, zoneIndex, coordIndexX, &min, &max, x_v.data());
    cgp_coord_write_data(file, baseIndex, zoneIndex, coordIndexY, &min, &max, y_v.data());
  } else {
    cgsize_t min = 0;
    cgsize_t max = 0;
    cgp_coord_write_data(file, baseIndex, zoneIndex, coordIndexX, &min, &max, NULL);
    cgp_coord_write_data(file, baseIndex, zoneIndex, coordIndexY, &min, &max, NULL);
  }

  // Write elements
  int sectionIndex;
  int start = 1;
  int end = sizes[1];

  cgp_section_write(file, baseIndex, zoneIndex, "GridElements", CGNS_ENUMV(TRI_3),
                    start, end, 0, &sectionIndex);

  if(rank == 0) {
    cgp_elements_write_data(file, baseIndex, zoneIndex, sectionIndex, start, end, cells.data());
  } else {
    cgp_elements_write_data(file, baseIndex, zoneIndex, sectionIndex, 0, 0, NULL);
  }

  int flowIndex;
  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, "FlowSolution0", CGNS_ENUMV(CellCenter), &flowIndex);

  cgsize_t numCells = mesh->cells->size * 9;
  cgsize_t minCell  = 9 * get_global_start_index(mesh->cells) + 1;
  cgsize_t maxCell  = minCell + numCells - 1;

  op_par_loop(save_values, "save_values", mesh->cells,
              op_arg_dat(data->Q[0][0],   -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

  int velXIndex;
  double *velX_data = getOP2Array(data->save_temp);
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", &velXIndex);
  cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, velXIndex, &minCell, &maxCell, velX_data);
  free(velX_data);

  op_par_loop(save_values, "save_values", mesh->cells,
              op_arg_dat(data->Q[0][1],   -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

  int velYIndex;
  double *velY_data = getOP2Array(data->save_temp);
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", &velYIndex);
  cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, velYIndex, &minCell, &maxCell, velY_data);
  free(velY_data);

  op_par_loop(save_values, "save_values", mesh->cells,
              op_arg_dat(data->p,         -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

  int pIndex;
  double *p_data = getOP2Array(data->save_temp);
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", &pIndex);
  cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, pIndex, &minCell, &maxCell, p_data);
  free(p_data);

  op_par_loop(save_values, "save_values", mesh->cells,
              op_arg_dat(data->vorticity, -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

  int vortIndex;
  double *vort_data = getOP2Array(data->save_temp);
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VorticityMagnitude", &vortIndex);
  cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, vortIndex, &minCell, &maxCell, vort_data);
  free(vort_data);

  int sIndex;
  if(ls) {
    op_par_loop(save_values, "save_values", mesh->cells,
                op_arg_dat(ls->step_s,      -1, OP_ID, 10, "double", OP_READ),
                op_arg_dat(data->save_temp, -1, OP_ID, 9, "double", OP_WRITE));

    double *s_data = getOP2Array(data->save_temp);
    cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Surface", &sIndex);
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, sIndex, &minCell, &maxCell, vort_data);
    free(s_data);
  }

  float exp[5];
  cg_goto(file, baseIndex, "end");
  cg_dataclass_write(CGNS_ENUMV(Dimensional));
  cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Degree));

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velXIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velYIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 1.0f; exp[1] = -1.0f; exp[2] = -2.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", pIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 0.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", vortIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  if(ls) {
    exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = 0.0f; exp[3] = 0.0f; exp[4] = 0.0f;
    cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", sIndex, "end");
    cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
  }

  cgp_close(file);

  if(rank == 0) {
    free(x_g);
    free(y_g);
  }
}

void save_solution_finalise(std::string filename, int numIter, double dt) {
  int rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  vector<double> times;
  char flowPtrs[numIter][32];

  for(int i = 0; i < numIter; i++) {
    times.push_back(i * dt);
    string name = "FlowSolution" + to_string(i);
    strcpy(flowPtrs[i], name.c_str());
  }

  int file;
  if (cgp_open(filename.c_str(), CG_MODE_MODIFY, &file)) {
    cgp_error_exit();
  }
  int baseIndex = 1;
  int zoneIndex = 1;

  // Create base iteration node
  cg_biter_write(file, baseIndex, "BaseIter", numIter);
  // Store time values of each iteration
  cg_gopath(file, "/Base/BaseIter");
  cgsize_t timeDims[1] = {times.size()};
  int timeIndex;
  cgp_array_write("TimeValues", CGNS_ENUMV(RealDouble), 1, timeDims, &timeIndex);
  if(rank == 0) {
    cgsize_t min = 1;
    cgsize_t max = times.size();
    cgp_array_write_data(timeIndex, &min, &max, times.data());
  } else {
    cgsize_t min = 0;
    cgsize_t max = 0;
    cgp_array_write_data(timeIndex, &min, &max, NULL);
  }

  // Create zone iteration node
  cg_ziter_write(file, baseIndex, zoneIndex, "ZoneIter");
  cg_gopath(file, "/Base/Zone/ZoneIter");
  cgsize_t flowPtrsDim[2] = {32, numIter};
  int flowPtrIndex;
  cgp_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, flowPtrsDim, &flowPtrIndex);
  if(rank == 0) {
    cgsize_t min[] = {1, 1};
    cgp_array_write_data(flowPtrIndex, min, flowPtrsDim, flowPtrs);
  } else {
    cgsize_t min[] = {0, 0};
    cgsize_t max[] = {0, 0};
    cgp_array_write_data(flowPtrIndex, min, max, NULL);
  }

  cg_simulation_type_write(file, baseIndex, CGNS_ENUMV(TimeAccurate));

  cgp_close(file);
}


void save_solution(std::string filename, DGMesh *mesh, INSData *data, int ind,
                   LS *ls, double finalTime, double nu) {
  // Calculate vorticity
  curl(mesh, data->Q[ind][0], data->Q[ind][1], data->vorticity);

  // Get local data (in same format that was originally passed to OP2)
  double *Ux   = (double *)malloc(10 * mesh->numCells * sizeof(double));
  double *Uy   = (double *)malloc(10 * mesh->numCells * sizeof(double));
  double *pr   = (double *)malloc(10 * mesh->numCells * sizeof(double));
  double *vort = (double *)malloc(10 * mesh->numCells * sizeof(double));
  double *x    = (double *)malloc(10 * mesh->numCells * sizeof(double));
  double *y    = (double *)malloc(10 * mesh->numCells * sizeof(double));
  double *s    = (double *)calloc(10 * mesh->numCells, sizeof(double));

  op_fetch_data(data->Q[ind][0], Ux);
  op_fetch_data(data->Q[ind][1], Uy);
  op_fetch_data(data->p, pr);
  op_fetch_data(data->vorticity, vort);
  op_fetch_data(mesh->x, x);
  op_fetch_data(mesh->y, y);
  if(ls) {
    op_fetch_data(ls->step_s, s);
  }

  // Gather onto root process
  int rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  double *Ux_g   = (double *)malloc(10 * mesh->numCells_g * sizeof(double));
  double *Uy_g   = (double *)malloc(10 * mesh->numCells_g * sizeof(double));
  double *pr_g   = (double *)malloc(10 * mesh->numCells_g * sizeof(double));
  double *vort_g = (double *)malloc(10 * mesh->numCells_g * sizeof(double));
  double *x_g    = (double *)malloc(10 * mesh->numCells_g * sizeof(double));
  double *y_g    = (double *)malloc(10 * mesh->numCells_g * sizeof(double));
  double *s_g    = (double *)malloc(10 * mesh->numCells_g * sizeof(double));

  gather_double_array(Ux_g, Ux, comm_size, mesh->numCells_g, mesh->numCells, 10);
  gather_double_array(Uy_g, Uy, comm_size, mesh->numCells_g, mesh->numCells, 10);
  gather_double_array(pr_g, pr, comm_size, mesh->numCells_g, mesh->numCells, 10);
  gather_double_array(vort_g, vort, comm_size, mesh->numCells_g, mesh->numCells, 10);
  gather_double_array(x_g, x, comm_size, mesh->numCells_g, mesh->numCells, 10);
  gather_double_array(y_g, y, comm_size, mesh->numCells_g, mesh->numCells, 10);
  gather_double_array(s_g, s, comm_size, mesh->numCells_g, mesh->numCells, 10);

  int file;
  if (cgp_open(filename.c_str(), CG_MODE_WRITE, &file)) {
    cgp_error_exit();
  }

  vector<double> x_v;
  vector<double> y_v;
  vector<double> u_v;
  vector<double> v_v;
  vector<double> pr_v;
  vector<double> vort_v;
  vector<double> s_v;
  vector<cgsize_t> cells(3 * mesh->numCells_g * 9);

  if(rank == 0) {
    if(DG_ORDER == 4) {
      get_data_vectors_order_4(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells,
                               Ux_g, Uy_g, pr_g, vort_g, x_g, y_g, s_g,
                               mesh->numCells_g);
    } else if(DG_ORDER == 3) {
      get_data_vectors_order_3(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells,
                               Ux_g, Uy_g, pr_g, vort_g, x_g, y_g, s_g,
                               mesh->numCells_g);
    } else if(DG_ORDER == 2) {
      get_data_vectors_order_2(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells,
                               Ux_g, Uy_g, pr_g, vort_g, x_g, y_g, s_g,
                               mesh->numCells_g);
    } else {
      get_data_vectors_order_1(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells,
                               Ux_g, Uy_g, pr_g, vort_g, x_g, y_g, s_g,
                               mesh->numCells_g);
    }
  }

  int baseIndex;
  int zoneIndex;
  string baseName = "Base";
  int cellDim = 2;
  int physicalDim = 2;
  cg_base_write(file, baseName.c_str(), cellDim, physicalDim, &baseIndex);
  // Create zone
  string zoneName = "Zone1";
  int sizes_i[3];
  if(rank == 0) {
    sizes_i[0] = x_v.size();
    sizes_i[1] = mesh->numCells_g * 9;
    sizes_i[2] = 0;
  }
  MPI_Bcast(sizes_i, 3, MPI_INT, 0, MPI_COMM_WORLD);
  cgsize_t sizes[] = {sizes_i[0], sizes_i[1], sizes_i[2]};
  cg_zone_write(file, baseIndex, zoneName.c_str(), sizes,
                CGNS_ENUMV(Unstructured), &zoneIndex);
  // Write grid coordinates
  int coordIndexX;
  int coordIndexY;
  cgp_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble), "CoordinateX", &coordIndexX);
  cgp_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble), "CoordinateY", &coordIndexY);

  if(rank == 0) {
    cgsize_t min = 1;
    cgsize_t max = x_v.size();
    cgp_coord_write_data(file, baseIndex, zoneIndex, coordIndexX, &min, &max, x_v.data());
    cgp_coord_write_data(file, baseIndex, zoneIndex, coordIndexY, &min, &max, y_v.data());
  } else {
    cgsize_t min = 0;
    cgsize_t max = 0;
    cgp_coord_write_data(file, baseIndex, zoneIndex, coordIndexX, &min, &max, NULL);
    cgp_coord_write_data(file, baseIndex, zoneIndex, coordIndexY, &min, &max, NULL);
  }

  // Write elements
  int sectionIndex;
  int start = 1;
  int end = sizes[1];

  cgp_section_write(file, baseIndex, zoneIndex, "GridElements", CGNS_ENUMV(TRI_3),
                    start, end, 0, &sectionIndex);

  if(rank == 0) {
    cgp_elements_write_data(file, baseIndex, zoneIndex, sectionIndex, start, end, cells.data());
  } else {
    cgp_elements_write_data(file, baseIndex, zoneIndex, sectionIndex, 0, 0, NULL);
  }

  int flowIndex;
  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, "FlowSolution", CGNS_ENUMV(Vertex), &flowIndex);

  // Write velocity x
  int velXIndex;
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", &velXIndex);
  if(rank == 0) {
    cgsize_t min = 1;
    cgsize_t max = u_v.size();
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, velXIndex, &min, &max, u_v.data());
  } else {
    cgsize_t min = 0;
    cgsize_t max = 0;
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, velXIndex, &min, &max, NULL);
  }

  // Write velocity y
  int velYIndex;
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", &velYIndex);
  if(rank == 0) {
    cgsize_t min = 1;
    cgsize_t max = v_v.size();
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, velYIndex, &min, &max, v_v.data());
  } else {
    cgsize_t min = 0;
    cgsize_t max = 0;
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, velYIndex, &min, &max, NULL);
  }

  // Write pressure
  int pIndex;
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", &pIndex);
  if(rank == 0) {
    cgsize_t min = 1;
    cgsize_t max = pr_v.size();
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, pIndex, &min, &max, pr_v.data());
  } else {
    cgsize_t min = 0;
    cgsize_t max = 0;
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, pIndex, &min, &max, NULL);
  }

  // Write vorticity
  int vortIndex;
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VorticityMagnitude", &vortIndex);
  if(rank == 0) {
    cgsize_t min = 1;
    cgsize_t max = vort_v.size();
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, vortIndex, &min, &max, vort_v.data());
  } else {
    cgsize_t min = 0;
    cgsize_t max = 0;
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, vortIndex, &min, &max, NULL);
  }

  // Write surface
  int sIndex;
  cgp_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Surface", &sIndex);
  if(rank == 0) {
    cgsize_t min = 1;
    cgsize_t max = s_v.size();
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, sIndex, &min, &max, s_v.data());
  } else {
    cgsize_t min = 0;
    cgsize_t max = 0;
    cgp_field_write_data(file, baseIndex, zoneIndex, flowIndex, sIndex, &min, &max, NULL);
  }

  float exp[5];
  cg_goto(file, baseIndex, "end");
  cg_dataclass_write(CGNS_ENUMV(Dimensional));
  cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Degree));

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velXIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velYIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 1.0f; exp[1] = -1.0f; exp[2] = -2.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", pIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 0.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", vortIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = 0.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", sIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  cgsize_t dim[2] = {1, 2};
  double infoData[] = {finalTime, nu};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("info");
  cg_gopath(file, "/Base/Zone1/info");
  int arrayIndex;
  cgp_array_write("info", CGNS_ENUMV(RealDouble), 2, dim, &arrayIndex);
  if(rank == 0) {
    cgsize_t min[] = {1, 1};
    cgsize_t max[] = {1, 2};
    cgp_array_write_data(arrayIndex, min, max, infoData);
  } else {
    cgsize_t min[] = {0, 0};
    cgsize_t max[] = {0, 0};
    cgp_array_write_data(arrayIndex, min, max, NULL);
  }

  cgp_close(file);

  free(s);
  free(Ux);
  free(Uy);
  free(pr);
  free(vort);
  free(x);
  free(y);
  free(s_g);
  free(Ux_g);
  free(Uy_g);
  free(pr_g);
  free(vort_g);
  free(x_g);
  free(y_g);
}