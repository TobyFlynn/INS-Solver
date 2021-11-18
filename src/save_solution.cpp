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
#include <algorithm>
#include <map>
#include <cmath>

#include "ins_data.h"
#include "dg_operators.h"

using namespace std;

struct Point {
  double x;
  double y;
  double u;
  double v;
  double pr;
  double vort;
  vector<int> cells;
  vector<int> pointNum;
  int counter;
};

struct cmpCoords {
    bool operator()(const pair<double,double>& a, const pair<double,double>& b) const {
        bool xCmp = abs(a.first - b.first) < 1e-8;
        bool yCmp = abs(a.second - b.second) < 1e-8;
        if(xCmp && yCmp) {
          return false;
        }
        return a < b;
    }
};

void get_data_vectors(INSData *data, int ind, vector<double> &x_v, vector<double> &y_v,
                      vector<double> &u_v, vector<double> &v_v, vector<double> &pr_v,
                      vector<double> &vort_v, vector<cgsize_t> &cells) {
  // Calculate vorticity
  curl(data->mesh, data->Q[ind][0], data->Q[ind][1], data->vorticity);

  // Get Data from OP2
  int numCells = op_get_size(data->mesh->cells);
  double *Ux   = (double *)malloc(15 * numCells * sizeof(double));
  double *Uy   = (double *)malloc(15 * numCells * sizeof(double));
  double *pr   = (double *)malloc(15 * numCells * sizeof(double));
  double *vort = (double *)malloc(15 * numCells * sizeof(double));
  double *x    = (double *)malloc(15 * numCells * sizeof(double));
  double *y    = (double *)malloc(15 * numCells * sizeof(double));

  op_fetch_data(data->Q[ind][0], Ux);
  op_fetch_data(data->Q[ind][1], Uy);
  op_fetch_data(data->p, pr);
  op_fetch_data(data->vorticity, vort);
  op_fetch_data(data->mesh->x, x);
  op_fetch_data(data->mesh->y, y);
  // Maps points to sub elements that they are part of.
  // Each line is 6 long (as 6 is the max number of sub elements within an original element that a point can be part of)
  // -1 is just padding to get each line to 6
  int cellMask[15][6] = {
    {0, -1, -1, -1, -1, -1}, // Point 0 is part of sub element 0
    {0, 1, 2, -1, -1, -1},    // 1
    {2, 3, 4, -1, -1, -1},    // 2
    {4, 5, 6, -1, -1, -1},    // 3
    {6, -1, -1, -1, -1, -1}, // End of first point row
    {0, 1, 7, -1, -1, -1},    // 5
    {1, 2, 3, 9, 8, 7},       // 6
    {3, 4, 5, 11, 10, 9},     // 7
    {6, 5, 11, -1, -1, -1}, // End of second point row
    {7, 8, 12, -1, -1, -1},   // 9
    {8, 9, 10, 14, 13, 12},   // 10
    {11, 10, 14, -1, -1, -1}, // End of third point row
    {12, 13, 15, -1, -1, -1}, // 12
    {14, 13, 15, -1, -1, -1}, // End of fourth point row
    {15, -1, -1, -1, -1, -1}  // 14
  };

  int pointNum[15][6] = {
    {0, -1, -1, -1, -1, -1}, // Point 0 is the first point of sub element 0
    {1, 0, 0, -1, -1, -1},  // 1
    {1, 0, 0, -1, -1, -1},  // 2
    {1, 0, 0, -1, -1, -1},  // 3
    {1, -1, -1, -1, -1, -1}, // End of first point row
    {2, 2, 0, -1, -1, -1},  // 5
    {1, 2, 2, 0, 0, 1},     // 6
    {1, 2, 2, 0, 0, 1},     // 7
    {2, 1, 1, -1, -1, -1}, // End of second point row
    {2, 2, 0, -1, -1, -1},  // 9
    {1, 2, 2, 0, 0, 1},     // 10
    {2, 1, 1, -1, -1, -1}, // End of third point row
    {2, 2, 0, -1, -1, -1},  // 12
    {2, 1, 1, -1, -1, -1}, // End of fourth point row
    {2, -1, -1, -1, -1, -1} // 14
  };

  // 16 sub elements per original element
  map<pair<double,double>,unique_ptr<Point>, cmpCoords> pointMap;
  for(int c = 0; c < numCells; c++) {
    int ind = c * 15;
    for(int p = 0; p < 15; p++) {
      pair<double,double> coords = make_pair(x[ind + p], y[ind + p]);
      unique_ptr<Point> point = make_unique<Point>();
      auto res = pointMap.insert(make_pair(coords, move(point)));
      if(res.second) {
        res.first->second->x    = x[ind + p];
        res.first->second->y    = y[ind + p];
        res.first->second->u    = Ux[ind + p];
        res.first->second->v    = Uy[ind + p];
        res.first->second->pr   = pr[ind + p];
        res.first->second->vort = vort[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 16 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter = 1;
      } else {
        res.first->second->u    += Ux[ind + p];
        res.first->second->v    += Uy[ind + p];
        res.first->second->pr   += pr[ind + p];
        res.first->second->vort += vort[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 16 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter++;
      }
    }
  }

  int index = 0;

  for(auto const &p : pointMap) {
    x_v.push_back(p.second->x);
    y_v.push_back(p.second->y);
    u_v.push_back(p.second->u / p.second->counter);
    v_v.push_back(p.second->v / p.second->counter);
    pr_v.push_back(p.second->pr / p.second->counter);
    vort_v.push_back(p.second->vort / p.second->counter);
    for(int i = 0; i < p.second->cells.size(); i++) {
      cells[p.second->cells[i] * 3 + p.second->pointNum[i]] = index + 1;
    }
    index++;
  }

  free(Ux);
  free(Uy);
  free(pr);
  free(vort);
  free(x);
  free(y);
}

void save_solution_iter(std::string filename, INSData *data, int ind, int iter) {
  int numCells = op_get_size(data->mesh->cells);
  vector<double> x_v;
  vector<double> y_v;
  vector<double> u_v;
  vector<double> v_v;
  vector<double> pr_v;
  vector<double> vort_v;
  vector<cgsize_t> cells(3 * numCells * 16);

  get_data_vectors(data, ind, x_v, y_v, u_v, v_v, pr_v, vort_v, cells);

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

  cg_close(file);
}

void save_solution_init(std::string filename, INSData *data) {
  int numCells = op_get_size(data->mesh->cells);
  vector<double> x_v;
  vector<double> y_v;
  vector<double> u_v;
  vector<double> v_v;
  vector<double> pr_v;
  vector<double> vort_v;
  vector<cgsize_t> cells(3 * numCells * 16);

  get_data_vectors(data, 0, x_v, y_v, u_v, v_v, pr_v, vort_v, cells);

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
  sizes[1] = numCells * 16;
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

  cg_close(file);
}

void save_solution_finalise(std::string filename, INSData *data, int numIter, double dt) {
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

void save_solution(std::string filename, INSData *data, int ind, double finalTime, double nu) {
  int numCells = op_get_size(data->mesh->cells);
  vector<double> x_v;
  vector<double> y_v;
  vector<double> u_v;
  vector<double> v_v;
  vector<double> pr_v;
  vector<double> vort_v;
  vector<cgsize_t> cells(3 * numCells * 16);

  get_data_vectors(data, ind, x_v, y_v, u_v, v_v, pr_v, vort_v, cells);

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
  sizes[1] = numCells * 16;
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

  cgsize_t dim[2] = {1, 2};
  double infoData[] = {finalTime, nu};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("info");
  cg_gopath(file, "/Base/Zone1/info");
  cg_array_write("info", CGNS_ENUMV(RealDouble), 2, dim, infoData);

  cg_close(file);
}
