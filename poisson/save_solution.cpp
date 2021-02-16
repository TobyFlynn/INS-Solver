// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

void save_solution(std::string filename, int numPts, int numCells, double *q0, int *cellMap) {
  vector<double> p(numPts);
  vector<double> num(numPts);

  int qNode0Ind = 0;
  int qNode1Ind = 4;
  int qNode2Ind = 14;

  for(int i = 0; i < numCells; i++) {
    int qCellInd = i * 15;
    int node0 = cellMap[i * 3];
    int node1 = cellMap[i * 3 + 1];
    int node2 = cellMap[i * 3 + 2];

    p[node0] += q0[qCellInd + qNode0Ind];
    num[node0]++;

    p[node1] += q0[qCellInd + qNode1Ind];
    num[node1]++;

    p[node2] += q0[qCellInd + qNode2Ind];
    num[node2]++;
  }

  for(int i = 0; i < numPts; i++) {
    p[i] = p[i] / num[i];
  }

  // Write out CGNS file
  int file;
  if (cg_open(filename.c_str(), CG_MODE_MODIFY, &file)) {
    cg_error_exit();
  }

  int baseIndex = 1;
  int zoneIndex = 1;
  int flowIndex;
  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, "FlowSolution", CGNS_ENUMV(Vertex), &flowIndex);

  // Write pressure
  int pIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", p.data(), &pIndex);

  float exp[5];
  cg_goto(file, baseIndex, "end");
  cg_dataclass_write(CGNS_ENUMV(Dimensional));
  cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Degree));

  exp[0] = 1.0f; exp[1] = -1.0f; exp[2] = -2.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", pIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
  cg_close(file);
}
