#include "load_mesh.h"

#include "pcgnslib.h"
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
  cgp_elements_read_data(file, baseIndex, zoneIndex, 1, start, end, cells.data());
  // CGNS starts numbering from 1 but OP2 starts from 0
  std::transform(cells.begin(), cells.end(), cgnsCells, [](int x) { return x - 1;});
}

template<>
void cgns_load_cells<long>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int start, int end) {
  std::vector<cgsize_t> cells((end - start + 1) * 3);
  cgsize_t parentData;
  cgp_elements_read_data(file, baseIndex, zoneIndex, 1, start, end, cells.data());
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
  int rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int *bcs_g = (int *)malloc(12 * sizeof(int));

  // Read CGNS grid
  int file;
  if(cgp_open(filename.c_str(), CG_MODE_READ, &file)) {
    cgp_error_exit();
  }

  int baseIndex = 1;
  int zoneIndex = 1;
  cgsize_t cg_numNodes;
  char zoneName[33];
  // Get zone name and size
  cg_zone_read(file, baseIndex, zoneIndex, zoneName, &cg_numNodes);
  *numNodes_g = (int) cg_numNodes;
  *numNodes   = compute_local_size(*numNodes_g, comm_size, rank);

  // Get vertices
  std::vector<double> x(*numNodes);
  std::vector<double> y(*numNodes);
  cgsize_t minVertex = compute_global_start(*numNodes_g, comm_size, rank) + 1;
  cgsize_t maxVertex = minVertex + *numNodes - 1;
  cgp_coord_read_data(file, baseIndex, zoneIndex, 1, &minVertex, &maxVertex, x.data());
  cgp_coord_read_data(file, baseIndex, zoneIndex, 2, &minVertex, &maxVertex, y.data());

  *coords_data = (double *)malloc(2 * (*numNodes) * sizeof(double));
  double *coords_ptr = *coords_data;
  for(int i = 0; i < x.size(); i++) {
    coords_ptr[2 * i] = x[i];
    coords_ptr[2 * i + 1] = y[i];
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
  *numCells_g = elementEnd - elementStart + 1;
  *numCells   = compute_local_size(*numCells_g, comm_size, rank);
  *cells_data = (int *)malloc((*numCells_g) * 3 * sizeof(int));
  int cellStart_g  = compute_global_start(*numCells_g, comm_size, rank) + 1;
  int cellEnd_g    = cellStart_g + (*numCells) - 1;
  cgns_load_cells<cgsize_t>(file, baseIndex, zoneIndex, *cells_data, cellStart_g, cellEnd_g);

  // Get edge data
  cg_gopath(file, "/Base/Zone1/Edges");
  char arrayName[33];
  DataType_t arrayDataType;
  int arrayRank;
  cgsize_t arrayDims[2];
  cg_array_info(1, arrayName, &arrayDataType, &arrayRank, arrayDims);
  *numEdges_g = arrayDims[1];
  *numEdges   = compute_local_size(*numEdges_g, comm_size, rank);
  *edge2node_data  = (int *)malloc(2 * (*numEdges) * sizeof(int));
  *edge2cell_data  = (int *)malloc(2 * (*numEdges) * sizeof(int));
  *edgeNum_data    = (int *)malloc(2 * (*numEdges) * sizeof(int));
  std::vector<int> edgeData(arrayDims[0] * (*numEdges));
  cgsize_t edgeStart[] = {1, compute_global_start(*numEdges_g, comm_size, rank) + 1};
  cgsize_t edgeEnd []  = {arrayDims[0], edgeStart[1] + (*numEdges) - 1};
  cgp_array_read_data(1, edgeStart, edgeEnd, edgeData.data());

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
  *numBoundaryEdges_g = barrayDims[1];
  *numBoundaryEdges   = compute_local_size(*numBoundaryEdges_g, comm_size, rank);
  *bedge2node_data = (int *)malloc(2 * (*numBoundaryEdges) * sizeof(int));
  *bedge2cell_data = (int *)malloc((*numBoundaryEdges) * sizeof(int));
  *bedgeNum_data   = (int *)malloc((*numBoundaryEdges) * sizeof(int));
  *bedge_type_data = (int *)malloc((*numBoundaryEdges) * sizeof(int));
  std::vector<int> bedgeData(barrayDims[0] * (*numBoundaryEdges));
  cgsize_t bedgeStart[] = {1, compute_global_start(*numBoundaryEdges_g, comm_size, rank) + 1};
  cgsize_t bedgeEnd []  = {barrayDims[0], bedgeStart[1] + (*numBoundaryEdges) - 1};
  cgp_array_read_data(1, bedgeStart, bedgeEnd, bedgeData.data());

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
  cgsize_t bcStart[] = {1, 1};
  cgp_array_read_data(1, bcStart, bcarrayDims, bcs_g);

  cg_close(file);

  for(int i = 0; i < 3; i++) {
    pressure_dirichlet[i] = bcs_g[0 + i];
    pressure_neumann[i] = bcs_g[3 + i];
    viscosity_dirichlet[i] = bcs_g[6 + i];
    viscosity_neumann[i] = bcs_g[9 + i];
  }

  free(bcs_g);
}
