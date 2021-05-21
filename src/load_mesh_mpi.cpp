#include "load_mesh.h"

#include "cgnslib.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include "ins_data.h"

#include "mpi.h"
#include "mpi_helper_func.h"

template<typename T>
void cgns_load_cells(int file, int baseIndex, int zoneIndex, int *cgnsCells, int start, int end);

template<>
void cgns_load_cells<int>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int start, int end) {
  std::vector<cgsize_t> cells((end - start + 1) * 3);
  cgsize_t parentData;
  cg_elements_partial_read(file, baseIndex, zoneIndex, 1, start, end, cells.data(), &parentData);
  // CGNS starts numbering from 1 but OP2 starts from 0
  std::transform(cells.begin(), cells.end(), cgnsCells, [](int x) { return x - 1;});
}

template<>
void cgns_load_cells<long>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int start, int end) {
  std::vector<cgsize_t> cells((end - start + 1) * 3);
  cgsize_t parentData;
  cg_elements_partial_read(file, baseIndex, zoneIndex, 1, start, end, cells.data(), &parentData);
  // CGNS starts numbering from 1 but OP2 starts from 0
  std::transform(cells.begin(), cells.end(), cgnsCells, [](long x) { return (int)x - 1;});
}

void load_mesh(std::string filename, INSData *data) {
  int rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  
  int *bcs_g = (int *)malloc(12 * sizeof(int));

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
  data->numNodes_g = (int) cg_numNodes;
  data->numNodes   = compute_local_size(data->numNodes_g, comm_size, rank);

  // Get vertices
  std::vector<double> x(data->numNodes);
  std::vector<double> y(data->numNodes);
  cgsize_t minVertex = compute_global_start(data->numNodes_g, comm_size, rank) + 1;
  cgsize_t maxVertex = minVertex + data->numNodes - 1;
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateX",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, x.data());
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateY",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, y.data());

  data->coords = (double *)malloc(2 * data->numNodes * sizeof(double));
  for(int i = 0; i < x.size(); i++) {
    data->coords[2 * i] = x[i];
    data->coords[2 * i + 1] = y[i];
  }

  // Get number of sections
  int numSections;
  cg_nsections(file, baseIndex, zoneIndex, &numSections);

  // Get cell section
  char sectionName[33];
  CGNS_ENUMT(ElementType_t) elementType;
  cgsize_t elementStart, elementEnd;
  int elementNumBoundary, parentFlag;
  cg_section_read(file, baseIndex, zoneIndex, 1, sectionName, &elementType,
                  &elementStart, &elementEnd, &elementNumBoundary, &parentFlag);

  // Get cells
  data->numCells_g = elementEnd - elementStart + 1;
  data->numCells   = compute_local_size(data->numCells_g, comm_size, rank);
  data->cgnsCells  = (int *)malloc(data->numCells_g * 3 * sizeof(int));
  int cellStart_g  = compute_global_start(data->numCells_g, comm_size, rank) + 1;
  int cellEnd_g    = cellStart_g + data->numCells - 1;
  cgns_load_cells<cgsize_t>(file, baseIndex, zoneIndex, data->cgnsCells, cellStart_g, cellEnd_g);

  // Get edge data
  cg_gopath(file, "/Base/Zone1/Edges");
  char arrayName[33];
  DataType_t arrayDataType;
  int arrayRank;
  cgsize_t arrayDims[2];
  cg_array_info(1, arrayName, &arrayDataType, &arrayRank, arrayDims);
  data->numEdges_g = arrayDims[1];
  data->numEdges   = compute_local_size(data->numEdges_g, comm_size, rank);
  data->edge2node_data  = (int *)malloc(2 * data->numEdges * sizeof(int));
  data->edge2cell_data  = (int *)malloc(2 * data->numEdges * sizeof(int));
  data->edgeNum_data    = (int *)malloc(2 * data->numEdges * sizeof(int));
  std::vector<int> edgeData(arrayDims[0] * data->numEdges);
  cgsize_t edgeStart[] = {1, compute_global_start(data->numEdges_g, comm_size, rank) + 1};
  cgsize_t edgeEnd []  = {arrayDims[0], edgeStart[1] + data->numEdges - 1};
  cgsize_t memDim   = arrayDims[0] * data->numEdges;
  cgsize_t memStart = 1;
  cgsize_t memEnd   = arrayDims[0] * data->numEdges;
  cg_array_general_read(1, edgeStart, edgeEnd, CGNS_ENUMV(Integer), 1, &memDim, &memStart, &memEnd, edgeData.data());

  for(int i = 0; i < data->numEdges; i++) {
    // - 1 as CGNS counts points from 1 but OP2 counts from 0
    // Cell index do start from one in this data
    data->edge2node_data[i * 2]     = edgeData[i * 6] - 1;
    data->edge2node_data[i * 2 + 1] = edgeData[i * 6 + 1] - 1;
    data->edge2cell_data[i * 2]     = edgeData[i * 6 + 2];
    data->edge2cell_data[i * 2 + 1] = edgeData[i * 6 + 3];
    data->edgeNum_data[i * 2]       = edgeData[i * 6 + 4];
    data->edgeNum_data[i * 2 + 1]   = edgeData[i * 6 + 5];
  }

  // Get boundary edge data
  cg_gopath(file, "/Base/Zone1/BoundaryEdges");
  char barrayName[33];
  DataType_t barrayDataType;
  int barrayRank;
  cgsize_t barrayDims[2];
  cg_array_info(1, barrayName, &barrayDataType, &barrayRank, barrayDims);
  data->numBoundaryEdges_g = barrayDims[1];
  data->numBoundaryEdges   = compute_local_size(data->numBoundaryEdges_g, comm_size, rank);
  data->bedge2node_data = (int *)malloc(2 * data->numBoundaryEdges * sizeof(int));
  data->bedge2cell_data = (int *)malloc(data->numBoundaryEdges * sizeof(int));
  data->bedgeNum_data   = (int *)malloc(data->numBoundaryEdges * sizeof(int));
  data->bedge_type_data = (int *)malloc(data->numBoundaryEdges * sizeof(int));
  std::vector<int> bedgeData(barrayDims[0] * data->numBoundaryEdges);
  cgsize_t bedgeStart[] = {1, compute_global_start(data->numBoundaryEdges_g, comm_size, rank) + 1};
  cgsize_t bedgeEnd []  = {barrayDims[0], bedgeStart[1] + data->numBoundaryEdges - 1};
  cgsize_t bmemDim   = barrayDims[0] * data->numBoundaryEdges;
  cgsize_t bmemStart = 1;
  cgsize_t bmemEnd   = barrayDims[0] * data->numBoundaryEdges;
  cg_array_general_read(1, bedgeStart, bedgeEnd, CGNS_ENUMV(Integer), 1, &bmemDim, &bmemStart, &bmemEnd, bedgeData.data());

  for(int i = 0; i < data->numBoundaryEdges; i++) {
    // - 1 as CGNS counts points from 1 but OP2 counts from 0
    // Cell index do start from one in this data
    data->bedge2node_data[i * 2]     = bedgeData[i * 5] - 1;
    data->bedge2node_data[i * 2 + 1] = bedgeData[i * 5 + 1] - 1;
    data->bedge2cell_data[i]         = bedgeData[i * 5 + 2];
    data->bedgeNum_data[i]           = bedgeData[i * 5 + 3];
    data->bedge_type_data[i]         = bedgeData[i * 5 + 4];
    if(data->bedge_type_data[i] < 0)
      std::cerr << "Error reading in boundary edge type" << std::endl;
  }

  cg_gopath(file, "/Base/Zone1/BCs");
  char bcarrayName[33];
  DataType_t bcarrayDataType;
  int bcarrayRank;
  cgsize_t bcarrayDims[2];
  cg_array_info(1, bcarrayName, &bcarrayDataType, &bcarrayRank, bcarrayDims);
  cg_array_read(1, bcs_g);

  cg_close(file);

  for(int i = 0; i < 3; i++) {
    data->pressure_dirichlet[i] = bcs_g[0 + i];
    data->pressure_neumann[i] = bcs_g[3 + i];
    data->viscosity_dirichlet[i] = bcs_g[6 + i];
    data->viscosity_neumann[i] = bcs_g[9 + i];
  }

  free(bcs_g);
}
