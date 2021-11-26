// Include OP2 stuff
#include "op_seq.h"
// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "dg_mesh.h"
#include "ins_data.h"
#include "ls.h"
#include "dg_operators.h"
#include "shared_save_functions.h"

using namespace std;

void save_solution_iter(std::string filename, DGMesh *mesh, INSData *data, int ind, LS *ls, int iter) {
  op_par_loop(save_order, "save_order", mesh->cells,
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);
  dats_to_update.push_back(ls->s);
  mesh->update_order(data->new_order, dats_to_update);

  int numCells = op_get_size(mesh->cells);
  vector<double> x_v;
  vector<double> y_v;
  vector<double> u_v;
  vector<double> v_v;
  vector<double> pr_v;
  vector<double> vort_v;
  vector<double> s_v;
  vector<cgsize_t> cells(3 * numCells * DG_SUB_CELLS);

  // Calculate vorticity
  curl(mesh, data->Q[ind][0], data->Q[ind][1], data->vorticity);

  // Get Data from OP2
  double *Ux_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *Uy_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *pr_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *vort_g = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *x_g    = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *y_g    = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *s_g    = (double *)calloc(DG_NP * numCells, sizeof(double));

  op_fetch_data(data->Q[ind][0], Ux_g);
  op_fetch_data(data->Q[ind][1], Uy_g);
  op_fetch_data(data->p, pr_g);
  op_fetch_data(data->vorticity, vort_g);
  op_fetch_data(mesh->x, x_g);
  op_fetch_data(mesh->y, y_g);
  if(ls) {
    op_fetch_data(ls->s, s_g);
  }

  if(DG_ORDER == 4) {
    get_data_vectors_order_4(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else if(DG_ORDER == 3) {
    get_data_vectors_order_3(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else if(DG_ORDER == 2) {
    get_data_vectors_order_2(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else {
    get_data_vectors_order_1(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  }

  free(Ux_g);
  free(Uy_g);
  free(pr_g);
  free(vort_g);
  free(x_g);
  free(y_g);
  free(s_g);

  int file;
  if (cg_open(filename.c_str(), CG_MODE_MODIFY, &file)) {
    cg_error_exit();
  }
  int baseIndex = 1;
  int zoneIndex = 1;

  int flowIndex;
  string flowName = "FlowSolution" + to_string(iter);
  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, flowName.c_str(), CGNS_ENUMV(Vertex), &flowIndex);

  // Write velocity x
  int velXIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", u_v.data(), &velXIndex);

  // Write velocity y
  int velYIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", v_v.data(), &velYIndex);

  // Write pressure
  int pIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", pr_v.data(), &pIndex);

  // Write vorticity
  int vortIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VorticityMagnitude", vort_v.data(), &vortIndex);

  // Write surface
  int sIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Surface", s_v.data(), &sIndex);

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

  cg_close(file);

  op_par_loop(ls_update_order, "ls_update_order", mesh->cells,
              op_arg_dat(mesh->order,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&ls->alpha,       1, "double", OP_READ),
              op_arg_dat(ls->s,           -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  mesh->update_order(data->new_order, dats_to_update);
}

void save_solution_init(std::string filename, DGMesh *mesh, INSData *data, LS *ls) {
  int numCells = op_get_size(mesh->cells);
  vector<double> x_v;
  vector<double> y_v;
  vector<double> u_v;
  vector<double> v_v;
  vector<double> pr_v;
  vector<double> vort_v;
  vector<double> s_v;
  vector<cgsize_t> cells(3 * numCells * DG_SUB_CELLS);

  // Calculate vorticity
  curl(mesh, data->Q[0][0], data->Q[0][1], data->vorticity);

  // Get Data from OP2
  double *Ux_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *Uy_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *pr_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *vort_g = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *x_g    = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *y_g    = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *s_g    = (double *)calloc(DG_NP * numCells, sizeof(double));

  op_fetch_data(data->Q[0][0], Ux_g);
  op_fetch_data(data->Q[0][1], Uy_g);
  op_fetch_data(data->p, pr_g);
  op_fetch_data(data->vorticity, vort_g);
  op_fetch_data(mesh->x, x_g);
  op_fetch_data(mesh->y, y_g);
  if(ls) {
    op_fetch_data(ls->s, s_g);
  }

  if(DG_ORDER == 4) {
    get_data_vectors_order_4(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else if(DG_ORDER == 3) {
    get_data_vectors_order_3(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else if(DG_ORDER == 2) {
    get_data_vectors_order_2(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else {
    get_data_vectors_order_1(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  }

  free(Ux_g);
  free(Uy_g);
  free(pr_g);
  free(vort_g);
  free(x_g);
  free(y_g);
  free(s_g);

  int file;
  if (cg_open(filename.c_str(), CG_MODE_WRITE, &file)) {
    cg_error_exit();
  }
  // Create base
  int baseIndex;
  int zoneIndex;
  int cellDim = 2;
  int physicalDim = 2;
  cg_base_write(file, "Base", cellDim, physicalDim, &baseIndex);

  // Create zone
  cgsize_t sizes[3];
  // Number of vertices
  sizes[0] = x_v.size();
  // Number of cells
  sizes[1] = numCells * DG_SUB_CELLS;
  // Number of boundary vertices (zero if elements not sorted)
  sizes[2] = 0;
  cg_zone_write(file, baseIndex, "Zone", sizes,
                CGNS_ENUMV(Unstructured), &zoneIndex);
  // Write grid coordinates
  int coordIndex;
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateX", x_v.data(), &coordIndex);
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateY", y_v.data(), &coordIndex);
  // Write elements
  int sectionIndex;
  int start = 1;
  int end = sizes[1];
  cg_section_write(file, baseIndex, zoneIndex, "GridElements", CGNS_ENUMV(TRI_3),
                   start, end, 0, cells.data(), &sectionIndex);

  // Write first flow solution
  int flowIndex;
  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, "FlowSolution0", CGNS_ENUMV(Vertex), &flowIndex);

  // Write velocity x
  int velXIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", u_v.data(), &velXIndex);

  // Write velocity y
  int velYIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", v_v.data(), &velYIndex);

  // Write pressure
  int pIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", pr_v.data(), &pIndex);

  // Write vorticity
  int vortIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VorticityMagnitude", vort_v.data(), &vortIndex);

  // Write surface
  int sIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Surface", s_v.data(), &sIndex);

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

  cg_close(file);
}

void save_solution_finalise(std::string filename, int numIter, double dt) {
  vector<double> times;
  char flowPtrs[numIter][32];

  for(int i = 0; i < numIter; i++) {
    times.push_back(i * dt);
    string name = "FlowSolution" + to_string(i);
    strcpy(flowPtrs[i], name.c_str());
  }

  int file;
  if (cg_open(filename.c_str(), CG_MODE_MODIFY, &file)) {
    cg_error_exit();
  }
  int baseIndex = 1;
  int zoneIndex = 1;

  // Create base iteration node
  cg_biter_write(file, baseIndex, "BaseIter", numIter);
  // Store time values of each iteration
  cg_gopath(file, "/Base/BaseIter");
  cgsize_t timeDims[1] = {times.size()};
  cg_array_write("TimeValues", CGNS_ENUMV(RealDouble), 1, timeDims, times.data());

  // Create zone iteration node
  cg_ziter_write(file, baseIndex, zoneIndex, "ZoneIter");
  cg_gopath(file, "/Base/Zone/ZoneIter");
  cgsize_t flowPtrsDim[2] = {32, numIter};
  cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, flowPtrsDim, flowPtrs);

  cg_simulation_type_write(file, baseIndex, CGNS_ENUMV(TimeAccurate));

  cg_close(file);
}

void save_solution(std::string filename, DGMesh *mesh, INSData *data, int ind, LS *ls, double finalTime, double nu) {
  op_par_loop(save_order, "save_order", mesh->cells,
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);
  dats_to_update.push_back(ls->s);
  mesh->update_order(data->new_order, dats_to_update);

  int numCells = op_get_size(mesh->cells);
  vector<double> x_v;
  vector<double> y_v;
  vector<double> u_v;
  vector<double> v_v;
  vector<double> pr_v;
  vector<double> vort_v;
  vector<double> s_v;
  vector<cgsize_t> cells(3 * numCells * DG_SUB_CELLS);

  // Calculate vorticity
  curl(mesh, data->Q[ind][0], data->Q[ind][1], data->vorticity);

  // Get Data from OP2
  double *Ux_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *Uy_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *pr_g   = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *vort_g = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *x_g    = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *y_g    = (double *)malloc(DG_NP * numCells * sizeof(double));
  double *s_g    = (double *)calloc(DG_NP * numCells, sizeof(double));

  op_fetch_data(data->Q[ind][0], Ux_g);
  op_fetch_data(data->Q[ind][1], Uy_g);
  op_fetch_data(data->p, pr_g);
  op_fetch_data(data->vorticity, vort_g);
  op_fetch_data(mesh->x, x_g);
  op_fetch_data(mesh->y, y_g);
  if(ls) {
    op_fetch_data(ls->s, s_g);
  }

  if(DG_ORDER == 4) {
    get_data_vectors_order_4(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else if(DG_ORDER == 3) {
    get_data_vectors_order_3(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else if(DG_ORDER == 2) {
    get_data_vectors_order_2(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  } else {
    get_data_vectors_order_1(x_v, y_v, u_v, v_v, pr_v, vort_v, s_v, cells, Ux_g,
                             Uy_g, pr_g, vort_g, x_g, y_g, s_g, numCells);
  }

  free(Ux_g);
  free(Uy_g);
  free(pr_g);
  free(vort_g);
  free(x_g);
  free(y_g);
  free(s_g);

  int file;
  if (cg_open(filename.c_str(), CG_MODE_WRITE, &file)) {
    cg_error_exit();
  }
  int baseIndex;
  int zoneIndex;
  string baseName = "Base";
  int cellDim = 2;
  int physicalDim = 2;
  cg_base_write(file, baseName.c_str(), cellDim, physicalDim, &baseIndex);
  // Create zone
  string zoneName = "Zone1";
  cgsize_t sizes[3];
  // Number of vertices
  sizes[0] = x_v.size();
  // Number of cells
  sizes[1] = numCells * DG_SUB_CELLS;
  // Number of boundary vertices (zero if elements not sorted)
  sizes[2] = 0;
  cg_zone_write(file, baseIndex, zoneName.c_str(), sizes,
                CGNS_ENUMV(Unstructured), &zoneIndex);
  // Write grid coordinates
  int coordIndex;
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateX", x_v.data(), &coordIndex);
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateY", y_v.data(), &coordIndex);

  // Write elements
  int sectionIndex;
  int start = 1;
  int end = sizes[1];

  cg_section_write(file, baseIndex, zoneIndex, "GridElements", CGNS_ENUMV(TRI_3),
                   start, end, 0, cells.data(), &sectionIndex);

  int flowIndex;
  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, "FlowSolution", CGNS_ENUMV(Vertex), &flowIndex);

  // Write velocity x
  int velXIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", u_v.data(), &velXIndex);

  // Write velocity y
  int velYIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", v_v.data(), &velYIndex);

  // Write pressure
  int pIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", pr_v.data(), &pIndex);

  // Write vorticity
  int vortIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VorticityMagnitude", vort_v.data(), &vortIndex);

  // Write surface
  int sIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Surface", s_v.data(), &sIndex);

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
  cg_array_write("info", CGNS_ENUMV(RealDouble), 2, dim, infoData);

  cg_close(file);

  op_par_loop(ls_update_order, "ls_update_order", mesh->cells,
              op_arg_dat(mesh->order,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&ls->alpha,       1, "double", OP_READ),
              op_arg_dat(ls->s,           -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  mesh->update_order(data->new_order, dats_to_update);
}
