#include "load_mesh.h"

#include "cgnslib.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include "ins_data.h"

template<typename T>
void cgns_load_cells(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells);

template<>
void cgns_load_cells<int>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells) {
  std::vector<cgsize_t> cells(numCells * 3);
  cgsize_t parentData;
  cg_elements_read(file, baseIndex, zoneIndex, 1, cells.data(), &parentData);
  // CGNS starts numbering from 1 but OP2 starts from 0
  std::transform(cells.begin(), cells.end(), cgnsCells, [](int x) { return x - 1;});
}

template<>
void cgns_load_cells<long>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells) {
  std::vector<cgsize_t> cells(numCells * 3);
  cgsize_t parentData;
  cg_elements_read(file, baseIndex, zoneIndex, 1, cells.data(), &parentData);
  // CGNS starts numbering from 1 but OP2 starts from 0
  std::transform(cells.begin(), cells.end(), cgnsCells, [](long x) { return (int)x - 1;});
}

void load_mesh(std::string filename, INSData *nsData) {
  // Read CGNS grid
  int file;
  if(cg_open(filename.c_str(), CG_MODE_READ, &file)) {
    cg_error_exit();
  }

  int baseIndex = 1;
  int zoneIndex = 1;
  cgsize_t cg_numNodes;
  char zoneName[33];
  // Get zone name and size
  cg_zone_read(file, baseIndex, zoneIndex, zoneName, &cg_numNodes);
  nsData->numNodes = (int) cg_numNodes;
  std::cout << "Zone Name: " << zoneName << std::endl;
  std::cout << "Zone Size: " << nsData->numNodes << std::endl;

  // Get vertices
  std::vector<double> x(nsData->numNodes);
  std::vector<double> y(nsData->numNodes);
  cgsize_t minVertex = 1;
  cgsize_t maxVertex = cg_numNodes;
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateX",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, x.data());
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateY",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, y.data());

  nsData->coords = (double *)malloc(2 * nsData->numNodes * sizeof(double));
  for(int i = 0; i < x.size(); i++) {
    nsData->coords[2 * i] = x[i];
    nsData->coords[2 * i + 1] = y[i];
  }

  // Get number of sections
  int numSections;
  cg_nsections(file, baseIndex, zoneIndex, &numSections);
  std::cout << "Number of sections: " << numSections << std::endl << std::endl;

  // Get cell section
  char sectionName[33];
  CGNS_ENUMT(ElementType_t) elementType;
  cgsize_t elementStart, elementEnd;
  int elementNumBoundary, parentFlag;
  cg_section_read(file, baseIndex, zoneIndex, 1, sectionName, &elementType,
                  &elementStart, &elementEnd, &elementNumBoundary, &parentFlag);
  std::cout << "Section 1: " << sectionName << std::endl;
  std::cout << "Element Type: " << ElementTypeName[elementType] << std::endl;
  std::cout << "Start: " << elementStart << " End: " << elementEnd << std::endl;
  std::cout << "Element Number Boundary: " << elementNumBoundary;
  std::cout << " Parent Flag: " << parentFlag << std::endl;

  // Get cells
  nsData->numCells = elementEnd - elementStart + 1;
  nsData->cgnsCells = (int *)malloc(nsData->numCells * 3 * sizeof(int));
  cgns_load_cells<cgsize_t>(file, baseIndex, zoneIndex, nsData->cgnsCells, nsData->numCells);

  // Get edge data
  cg_gopath(file, "/Base/Zone1/Edges");
  char arrayName[33];
  DataType_t arrayDataType;
  int arrayRank;
  cgsize_t arrayDims[2];
  cg_array_info(1, arrayName, &arrayDataType, &arrayRank, arrayDims);
  std::cout << "Array Name: " << arrayName << std::endl;
  std::cout << "Array Dims: " << arrayDims[0] << " " << arrayDims[1] << std::endl;
  nsData->numEdges = arrayDims[1];
  nsData->edge2node_data = (int *)malloc(2 * nsData->numEdges * sizeof(int));
  nsData->edge2cell_data = (int *)malloc(2 * nsData->numEdges * sizeof(int));
  nsData->edgeNum_data   = (int *)malloc(2 * nsData->numEdges * sizeof(int));
  std::vector<int> edgeData(arrayDims[0] * arrayDims[1]);
  cg_array_read(1, edgeData.data());

  for(int i = 0; i < nsData->numEdges; i++) {
    // - 1 as CGNS counts points from 1 but OP2 counts from 0
    // Cell index do start from one in this data
    nsData->edge2node_data[i * 2]     = edgeData[i * 6] - 1;
    nsData->edge2node_data[i * 2 + 1] = edgeData[i * 6 + 1] - 1;
    nsData->edge2cell_data[i * 2]     = edgeData[i * 6 + 2];
    nsData->edge2cell_data[i * 2 + 1] = edgeData[i * 6 + 3];
    nsData->edgeNum_data[i * 2]       = edgeData[i * 6 + 4];
    nsData->edgeNum_data[i * 2 + 1]   = edgeData[i * 6 + 5];
  }

  // Get boundary edge data
  cg_gopath(file, "/Base/Zone1/BoundaryEdges");
  char barrayName[33];
  DataType_t barrayDataType;
  int barrayRank;
  cgsize_t barrayDims[2];
  cg_array_info(1, barrayName, &barrayDataType, &barrayRank, barrayDims);
  std::cout << "Array Name: " << barrayName << std::endl;
  std::cout << "Array Dims: " << barrayDims[0] << " " << barrayDims[1] << std::endl;
  nsData->numBoundaryEdges = barrayDims[1];
  nsData->bedge2node_data = (int *)malloc(2 * nsData->numBoundaryEdges * sizeof(int));
  nsData->bedge2cell_data = (int *)malloc(nsData->numBoundaryEdges * sizeof(int));
  nsData->bedgeNum_data   = (int *)malloc(nsData->numBoundaryEdges * sizeof(int));
  nsData->bedge_type_data = (int *)malloc(nsData->numBoundaryEdges * sizeof(int));
  std::vector<int> bedgeData(barrayDims[0] * barrayDims[1]);
  cg_array_read(1, bedgeData.data());

  for(int i = 0; i < nsData->numBoundaryEdges; i++) {
    // - 1 as CGNS counts points from 1 but OP2 counts from 0
    // Cell index do start from one in this data
    nsData->bedge2node_data[i * 2]     = bedgeData[i * 5] - 1;
    nsData->bedge2node_data[i * 2 + 1] = bedgeData[i * 5 + 1] - 1;
    nsData->bedge2cell_data[i]         = bedgeData[i * 5 + 2];
    nsData->bedgeNum_data[i]           = bedgeData[i * 5 + 3];
    nsData->bedge_type_data[i]         = bedgeData[i * 5 + 4];
    if(nsData->bedge_type_data[i] < 0)
      std::cerr << "Error reading in boundary edge type" << std::endl;
  }

  cg_close(file);
}
