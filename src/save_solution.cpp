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

void save_solution(std::string filename, int numPts, int numCells, double *q0,
                   double *q1, double *p, int *cellMap) {
  vector<double> velX(numPts);
  vector<double> velY(numPts);
  vector<double> pr(numPts);
  vector<int> num(numPts);

  for(int i = 0; i < numPts; i++) {
    velX[i] = 0.0;
    velY[i] = 0.0;
    pr[i] = 0.0;
    num[i] = 0;
  }

  int qNode0Ind = 0;
  int qNode1Ind = 4;
  int qNode2Ind = 14;

  for(int i = 0; i < numCells; i++) {
    int qCellInd = i * 15;
    int node0 = cellMap[i * 3];
    int node1 = cellMap[i * 3 + 1];
    int node2 = cellMap[i * 3 + 2];

    velX[node0] += q0[qCellInd + qNode0Ind];
    velY[node0] += q1[qCellInd + qNode0Ind];
    pr[node0] += p[qCellInd + qNode0Ind];
    num[node0]++;

    velX[node1] += q0[qCellInd + qNode1Ind];
    velY[node1] += q1[qCellInd + qNode1Ind];
    pr[node1] += p[qCellInd + qNode1Ind];
    num[node1]++;

    velX[node2] += q0[qCellInd + qNode2Ind];
    velY[node2] += q1[qCellInd + qNode2Ind];
    pr[node2] += p[qCellInd + qNode2Ind];
    num[node2]++;
  }

  for(int i = 0; i < numPts; i++) {
    velX[i] = velX[i] / (double)num[i];
    velY[i] = velY[i] / (double)num[i];
    pr[i] = pr[i] / (double)num[i];
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
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", pr.data(), &pIndex);

  // Write velocity x
  int velXIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", velX.data(), &velXIndex);

  // Write velocity y
  int velYIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", velY.data(), &velYIndex);

  float exp[5];
  cg_goto(file, baseIndex, "end");
  cg_dataclass_write(CGNS_ENUMV(Dimensional));
  cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Degree));

  exp[0] = 1.0f; exp[1] = -1.0f; exp[2] = -2.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", pIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velXIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velYIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  cg_close(file);
}

void save_solution_cell(std::string filename, int numPts, int numCells, double *q0,
                        double *q1, int *cellMap) {
  vector<double> velX(numCells);
  vector<double> velY(numCells);

  int qNode0Ind = 0;
  int qNode1Ind = 4;
  int qNode2Ind = 14;

  for(int i = 0; i < numCells; i++) {
    int qCellInd = i * 15;
    int node0 = cellMap[i * 3];
    int node1 = cellMap[i * 3 + 1];
    int node2 = cellMap[i * 3 + 2];

    velX[i] = 0.0;
    velY[i] = 0.0;

    for(int j = 0; j < 15; j++) {
      velX[i] += q0[qCellInd + j];
      velY[i] += q1[qCellInd + j];
    }

    velX[i] = velX[i] / 15.0;
    velY[i] = velY[i] / 15.0;
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
  cg_sol_write(file, baseIndex, zoneIndex, "FlowSolution", CGNS_ENUMV(CellCenter), &flowIndex);

  // Write velocity x
  int velXIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", velX.data(), &velXIndex);

  // Write velocity y
  int velYIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", velY.data(), &velYIndex);

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

void save_solution_t(std::string filename, int numPts, int numCells, double *q0, double *q1,
                     double *p, double *pRHS, double *px, double *py, double *utx, double *uty,
                     double *uttx, double *utty, double *visx, double *visy, int *cellMap) {
  vector<double> velX(numPts);
  vector<double> velY(numPts);
  vector<double> pr(numPts);
  vector<double> prRHS(numPts);
  vector<double> prX(numPts);
  vector<double> prY(numPts);
  vector<double> utX(numPts);
  vector<double> utY(numPts);
  vector<double> uttX(numPts);
  vector<double> uttY(numPts);
  vector<double> visX(numPts);
  vector<double> visY(numPts);
  vector<int> num(numPts);

  for(int i = 0; i < numPts; i++) {
    velX[i] = 0.0;
    velY[i] = 0.0;
    pr[i] = 0.0;
    prRHS[i] = 0.0;
    prX[i] = 0.0;
    prY[i] = 0.0;
    utX[i] = 0.0;
    utY[i] = 0.0;
    uttX[i] = 0.0;
    uttY[i] = 0.0;
    visX[i] = 0.0;
    visY[i] = 0.0;
    num[i] = 0;
  }

  int qNode0Ind = 0;
  int qNode1Ind = 4;
  int qNode2Ind = 14;

  for(int i = 0; i < numCells; i++) {
    int qCellInd = i * 15;
    int node0 = cellMap[i * 3];
    int node1 = cellMap[i * 3 + 1];
    int node2 = cellMap[i * 3 + 2];

    velX[node0] += q0[qCellInd + qNode0Ind];
    velY[node0] += q1[qCellInd + qNode0Ind];
    pr[node0] += p[qCellInd + qNode0Ind];
    prRHS[node0] += pRHS[qCellInd + qNode0Ind];
    prX[node0] += px[qCellInd + qNode0Ind];
    prY[node0] += py[qCellInd + qNode0Ind];
    utX[node0] += utx[qCellInd + qNode0Ind];
    utY[node0] += uty[qCellInd + qNode0Ind];
    uttX[node0] += uttx[qCellInd + qNode0Ind];
    uttY[node0] += utty[qCellInd + qNode0Ind];
    visX[node0] += visx[qCellInd + qNode0Ind];
    visY[node0] += visy[qCellInd + qNode0Ind];
    num[node0]++;

    velX[node1] += q0[qCellInd + qNode1Ind];
    velY[node1] += q1[qCellInd + qNode1Ind];
    pr[node1] += p[qCellInd + qNode1Ind];
    prRHS[node1] += pRHS[qCellInd + qNode1Ind];
    prX[node1] += px[qCellInd + qNode1Ind];
    prY[node1] += py[qCellInd + qNode1Ind];
    utX[node1] += utx[qCellInd + qNode1Ind];
    utY[node1] += uty[qCellInd + qNode1Ind];
    uttX[node1] += uttx[qCellInd + qNode1Ind];
    uttY[node1] += utty[qCellInd + qNode1Ind];
    visX[node1] += visx[qCellInd + qNode1Ind];
    visY[node1] += visy[qCellInd + qNode1Ind];
    num[node1]++;

    velX[node2] += q0[qCellInd + qNode2Ind];
    velY[node2] += q1[qCellInd + qNode2Ind];
    pr[node2] += p[qCellInd + qNode2Ind];
    prRHS[node2] += pRHS[qCellInd + qNode2Ind];
    prX[node2] += px[qCellInd + qNode2Ind];
    prY[node2] += py[qCellInd + qNode2Ind];
    utX[node2] += utx[qCellInd + qNode2Ind];
    utY[node2] += uty[qCellInd + qNode2Ind];
    uttX[node2] += uttx[qCellInd + qNode2Ind];
    uttY[node2] += utty[qCellInd + qNode2Ind];
    visX[node2] += visx[qCellInd + qNode2Ind];
    visY[node2] += visy[qCellInd + qNode2Ind];
    num[node2]++;
  }

  for(int i = 0; i < numPts; i++) {
    velX[i] = velX[i] / (double)num[i];
    velY[i] = velY[i] / (double)num[i];
    pr[i] = pr[i] / (double)num[i];
    prRHS[i] = prRHS[i] / (double)num[i];
    prX[i] = prX[i] / (double)num[i];
    prY[i] = prY[i] / (double)num[i];
    utX[i] = utX[i] / (double)num[i];
    utY[i] = utY[i] / (double)num[i];
    uttX[i] = uttX[i] / (double)num[i];
    uttY[i] = uttY[i] / (double)num[i];
    visX[i] = visX[i] / (double)num[i];
    visY[i] = visY[i] / (double)num[i];
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
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", pr.data(), &pIndex);

  // Write velocity x
  int velXIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", velX.data(), &velXIndex);

  // Write velocity y
  int velYIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", velY.data(), &velYIndex);

  int tIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "pRHS", prRHS.data(), &tIndex);
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "pX", prX.data(), &tIndex);
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "pY", prY.data(), &tIndex);
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "utX", utX.data(), &tIndex);
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "utY", utY.data(), &tIndex);
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "uttX", uttX.data(), &tIndex);
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "uttY", uttY.data(), &tIndex);
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "visX", visX.data(), &tIndex);
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "visY", visY.data(), &tIndex);

  cg_close(file);
}
