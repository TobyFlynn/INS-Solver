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

void load_mesh(std::string filename, double **coords_data, int **cells_data,
               int **edge2node_data, int **edge2cell_data,
               int **bedge2node_data, int **bedge2cell_data,
               int **bedge_type_data, int **edgeNum_data, int **bedgeNum_data,
               int *numNodes_g, int *numCells_g, int *numEdges_g,
               int *numBoundaryEdges_g, int *numNodes, int *numCells,
               int *numEdges, int *numBoundaryEdges, int *pressure_dirichlet,
               int *pressure_neumann, int *viscosity_dirichlet,
               int *viscosity_neumann) {
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
  *numNodes = (int) cg_numNodes;
  std::cout << "Zone Name: " << zoneName << std::endl;
  std::cout << "Zone Size: " << *numNodes << std::endl;

  // Get vertices
  std::vector<double> x(*numNodes);
  std::vector<double> y(*numNodes);
  cgsize_t minVertex = 1;
  cgsize_t maxVertex = cg_numNodes;
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateX",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, x.data());
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateY",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, y.data());

  *coords_data = (double *)malloc(2 * (*numNodes) * sizeof(double));
  double *coords_ptr = *coords_data;
  for(int i = 0; i < x.size(); i++) {
    coords_ptr[2 * i] = x[i];
    coords_ptr[2 * i + 1] = y[i];
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
  *numCells = elementEnd - elementStart + 1;
  *cells_data = (int *)malloc((*numCells) * 3 * sizeof(int));
  cgns_load_cells<cgsize_t>(file, baseIndex, zoneIndex, *cells_data, *numCells);

  // Get edge data
  cg_gopath(file, "/Base/Zone1/Edges");
  char arrayName[33];
  DataType_t arrayDataType;
  int arrayRank;
  cgsize_t arrayDims[2];
  cg_array_info(1, arrayName, &arrayDataType, &arrayRank, arrayDims);
  std::cout << "Array Name: " << arrayName << std::endl;
  std::cout << "Array Dims: " << arrayDims[0] << " " << arrayDims[1] << std::endl;
  *numEdges = arrayDims[1];
  *edge2node_data = (int *)malloc(2 * (*numEdges) * sizeof(int));
  *edge2cell_data = (int *)malloc(2 * (*numEdges) * sizeof(int));
  *edgeNum_data   = (int *)malloc(2 * (*numEdges) * sizeof(int));
  std::vector<int> edgeData(arrayDims[0] * arrayDims[1]);
  cg_array_read(1, edgeData.data());

  int *edge2node_ptr = *edge2node_data;
  int *edge2cell_ptr = *edge2cell_data;
  int *edgeNum_ptr   = *edgeNum_data;

  for(int i = 0; i < (*numEdges); i++) {
    // - 1 as CGNS counts points from 1 but OP2 counts from 0
    // Cell index do start from one in this data
    edge2node_ptr[i * 2]     = edgeData[i * 6] - 1;
    edge2node_ptr[i * 2 + 1] = edgeData[i * 6 + 1] - 1;
    edge2cell_ptr[i * 2]     = edgeData[i * 6 + 2];
    edge2cell_ptr[i * 2 + 1] = edgeData[i * 6 + 3];
    edgeNum_ptr[i * 2]       = edgeData[i * 6 + 4];
    edgeNum_ptr[i * 2 + 1]   = edgeData[i * 6 + 5];
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
  *numBoundaryEdges = barrayDims[1];
  *bedge2node_data = (int *)malloc(2 * (*numBoundaryEdges) * sizeof(int));
  *bedge2cell_data = (int *)malloc((*numBoundaryEdges) * sizeof(int));
  *bedgeNum_data   = (int *)malloc((*numBoundaryEdges) * sizeof(int));
  *bedge_type_data = (int *)malloc((*numBoundaryEdges) * sizeof(int));
  std::vector<int> bedgeData(barrayDims[0] * barrayDims[1]);
  cg_array_read(1, bedgeData.data());

  int *bedge2node_ptr = *bedge2node_data;
  int *bedge2cell_ptr = *bedge2cell_data;
  int *bedgeNum_ptr   = *bedgeNum_data;
  int *bedge_type_ptr = *bedge_type_data;

  for(int i = 0; i < (*numBoundaryEdges); i++) {
    // - 1 as CGNS counts points from 1 but OP2 counts from 0
    // Cell index do start from one in this data
    bedge2node_ptr[i * 2]     = bedgeData[i * 5] - 1;
    bedge2node_ptr[i * 2 + 1] = bedgeData[i * 5 + 1] - 1;
    bedge2cell_ptr[i]         = bedgeData[i * 5 + 2];
    bedgeNum_ptr[i]           = bedgeData[i * 5 + 3];
    bedge_type_ptr[i]         = bedgeData[i * 5 + 4];
    if(bedge_type_ptr[i] < 0)
      std::cerr << "Error reading in boundary edge type" << std::endl;
  }

  cg_gopath(file, "/Base/Zone1/BCs");
  char bcarrayName[33];
  DataType_t bcarrayDataType;
  int bcarrayRank;
  cgsize_t bcarrayDims[2];
  cg_array_info(1, bcarrayName, &bcarrayDataType, &bcarrayRank, bcarrayDims);
  std::cout << "Array Name: " << bcarrayName << std::endl;
  std::cout << "Array Dims: " << bcarrayDims[0] << " " << bcarrayDims[1] << std::endl;
  std::vector<int> bcData(bcarrayDims[0] * bcarrayDims[1]);
  cg_array_read(1, bcData.data());
  for(int i = 0; i < 3; i++) {
    pressure_dirichlet[i] = bcData[0 + i];
    pressure_neumann[i] = bcData[3 + i];
    viscosity_dirichlet[i] = bcData[6 + i];
    viscosity_neumann[i] = bcData[9 + i];
  }

  cg_close(file);

  for(int i = 0; i < x.size() * 2; i++) {
    if(isnan(coords_ptr[i]))
      printf("***NaN***\n");
  }
}
