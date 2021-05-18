// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <cmath>
#include <getopt.h>
#include <cstring>

using namespace std;

// Stuff for parsing command line arguments
extern char *optarg;
extern int  optind, opterr, optopt;
static struct option options[] = {
  {"file", required_argument, 0, 0},
  {0,    0,                  0,  0}
};

int main(int argc, char **argv) {
  string fileName = "end.cgns";
  int opt_index = 0;
  while(getopt_long_only(argc, argv, "", options, &opt_index) != -1) {
    if(strcmp((char*)options[opt_index].name,"file") == 0) fileName = optarg;
  }

  int file;
  if (cg_open(fileName.c_str(), CG_MODE_READ, &file)) {
    cg_error_exit();
  }

  int baseIndex = 1;
  int zoneIndex = 1;
  cgsize_t cg_numNodes;
  char zoneName[33];
  // Get zone name and size
  cg_zone_read(file, baseIndex, zoneIndex, zoneName, &cg_numNodes);
  int numNodes = (int) cg_numNodes;

  vector<double> x(numNodes);
  vector<double> y(numNodes);

  cgsize_t minVertex = 1;
  cgsize_t maxVertex = cg_numNodes;
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateX",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, x.data());
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateY",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, y.data());

  CGNS_ENUMT(ElementType_t) elementType;
  cgsize_t elementStart, elementEnd;
  int elementNumBoundary, parentFlag;
  char sectionName[33];
  cg_section_read(file, baseIndex, zoneIndex, 1, sectionName, &elementType,
                  &elementStart, &elementEnd, &elementNumBoundary, &parentFlag);
  // Get cells
  int numCells = elementEnd - elementStart + 1;
  vector<cgsize_t> cells(3 * numCells);
  cgsize_t parentData;
  cg_elements_read(file, baseIndex, zoneIndex, 1, cells.data(), &parentData);

  vector<double> u(numNodes);
  vector<double> v(numNodes);
  vector<double> pr(numNodes);
  int flowIndex = 1;
  cg_field_read(file, baseIndex, zoneIndex, flowIndex, "VelocityX", CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, u.data());
  cg_field_read(file, baseIndex, zoneIndex, flowIndex, "VelocityY", CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, v.data());
  cg_field_read(file, baseIndex, zoneIndex, flowIndex, "Pressure", CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, pr.data());

  cg_gopath(file, "/Base/Zone1/info");
  double data[2];
  cg_array_read(1, data);
  double time = data[0];
  double nu = data[1];

  cg_close(file);

  vector<double> uErr(numNodes);
  vector<double> vErr(numNodes);
  vector<double> prErr(numNodes);

  const double PI = 3.141592653589793238463;
  for(int i = 0; i < numNodes; i++) {
    double uExact  = -sin(2.0 * PI * y[i]) * exp(-nu * 4.0 * PI * PI * time);
    double vExact  = sin(2.0 * PI * x[i]) * exp(-nu * 4.0 * PI * PI * time);
    double prExact = -cos(2.0 * PI * x[i]) * cos(2.0 * PI * y[i]) * exp(-nu * 8.0 * PI * PI * time);

    uErr[i]  = abs(uExact - u[i]);
    vErr[i]  = abs(vExact - v[i]);
    prErr[i] = abs(prExact - pr[i]);
  }

  if (cg_open("err.cgns", CG_MODE_WRITE, &file)) {
    cg_error_exit();
  }

  string baseName = "Base";
  int cellDim = 2;
  int physicalDim = 2;
  cg_base_write(file, baseName.c_str(), cellDim, physicalDim, &baseIndex);
  // Create zone
  string zoneName2 = "Zone";
  cgsize_t sizes[3];
  // Number of vertices
  sizes[0] = x.size();
  // Number of cells
  sizes[1] = numCells;
  // Number of boundary vertices (zero if elements not sorted)
  sizes[2] = 0;
  cg_zone_write(file, baseIndex, zoneName2.c_str(), sizes,
                CGNS_ENUMV(Unstructured), &zoneIndex);
  // Write grid coordinates
  int coordIndex;
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateX", x.data(), &coordIndex);
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateY", y.data(), &coordIndex);

  // Write elements
  int sectionIndex;
  int start = 1;
  int end = sizes[1];

  cg_section_write(file, baseIndex, zoneIndex, "GridElements", CGNS_ENUMV(TRI_3),
                   start, end, 0, cells.data(), &sectionIndex);

  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, "FlowSolution", CGNS_ENUMV(Vertex), &flowIndex);

  // Write velocity x
  int velXIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", uErr.data(), &velXIndex);

  // Write velocity y
  int velYIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", vErr.data(), &velYIndex);

  // Write pressure
  int pIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", prErr.data(), &pIndex);

  cg_close(file);
}
