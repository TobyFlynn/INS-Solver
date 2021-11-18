#ifndef __INS_LOADMESH_H
#define __INS_LOADMESH_H

#include <string>

#include "ins_data.h"

void load_mesh(std::string filename, double **coords_data, int **cells_data,
               int **edge2node_data, int **edge2cell_data,
               int **bedge2node_data, int **bedge2cell_data,
               int **bedge_type_data, int **edgeNum_data, int **bedgeNum_data,
               int *numNodes_g, int *numCells_g, int *numEdges_g,
               int *numBoundaryEdges_g, int *numNodes, int *numCells,
               int *numEdges, int *numBoundaryEdges, int *pressure_dirichlet,
               int *pressure_neumann, int *viscosity_dirichlet,
               int *viscosity_neumann);

#endif
