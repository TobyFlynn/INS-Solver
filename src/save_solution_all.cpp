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

using namespace std;

struct Point {
  double x;
  double y;
  double u;
  double v;
  vector<int> cells;
  vector<int> pointNum;
};

struct cmpCoords {
    bool operator()(const pair<double,double>& a, const pair<double,double>& b) const {
        bool xCmp = abs(a.first - b.first) < 1e-12;
        bool yCmp = abs(a.second - b.second) < 1e-12;
        if(xCmp && yCmp) {
          return false;
        }
        return a < b;
    }
};

void save_solution_all(std::string filename, int numCells, double *q0, double *q1, double *x, double *y) {

  int mask[] = {0, 3, 4, 5, 1, 11, 12, 13, 6, 10, 14, 7, 9, 8, 2};

  map<pair<double,double>,unique_ptr<Point>, cmpCoords> pointMap;
  for(int i = 0; i < numCells; i++) {
    int ind = i * 15;
    for(int j = 0; j < 15; j++) {
      pair<double,double> coords = make_pair(x[ind + j], y[ind + j]);
      if(pointMap.count(coords) == 0) {
        unique_ptr<Point> p = make_unique<Point>();
        p->x = x[ind + j];
        p->y = y[ind + j];
        p->u = q0[ind + j];
        p->v = q1[ind + j];
        p->cells.push_back(i);
        p->pointNum.push_back(mask[j]);
        pointMap.insert(make_pair(coords, move(p)));
      } else {
        pointMap.at(coords)->u += q0[ind + j];
        pointMap.at(coords)->v += q1[ind + j];
        pointMap.at(coords)->cells.push_back(i);
        pointMap.at(coords)->pointNum.push_back(mask[j]);
      }
    }
  }

  vector<double> x_v;
  vector<double> y_v;
  vector<double> u_v;
  vector<double> v_v;
  vector<cgsize_t> cells(15 * numCells);
  int index = 0;

  for(auto const &p : pointMap) {
    x_v.push_back(p.second->x);
    y_v.push_back(p.second->y);
    u_v.push_back(p.second->u / cells.size());
    v_v.push_back(p.second->v / cells.size());
    for(int i = 0; i < p.second->cells.size(); i++) {
      cells[p.second->cells[i] * 15 + p.second->pointNum[i]] = index + 1;
    }
    index++;
  }

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
  string zoneName = "Zone";
  cgsize_t sizes[3];
  // Number of vertices
  sizes[0] = x_v.size();
  // Number of cells
  sizes[1] = numCells;
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
  cg_section_write(file, baseIndex, zoneIndex, "GridElements", CGNS_ENUMV(TRI_15),
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

  cg_close(file);
}
