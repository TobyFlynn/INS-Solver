#include "ins_data.h"

#include <string>
#include <memory>

using namespace std;

INSData::INSData() {

}

INSData::~INSData() {
  free(coords);
  free(cgnsCells);
  free(edge2node_data);
  free(edge2cell_data);
  free(bedge2node_data);
  free(bedge2cell_data);
  free(bedge_type_data);
  free(edgeNum_data);
  free(bedgeNum_data);

  free(nodeX_data);
  free(nodeY_data);
  free(x_data);
  free(y_data);
  free(xr_data);
  free(yr_data);
  free(xs_data);
  free(ys_data);
  free(rx_data);
  free(ry_data);
  free(sx_data);
  free(sy_data);
  free(nx_data);
  free(ny_data);
  free(fscale_data);
  for(int i = 0; i < 3; i++) {
    free(Q_data[i]);
  }
  for(int i = 0; i < 4; i++) {
    free(F_data[i]);
  }
}

void INSData::initOP2() {
  // Initialise memory
  nodeX_data = (double*)malloc(3 * numCells * sizeof(double));
  nodeY_data = (double*)malloc(3 * numCells * sizeof(double));
  x_data  = (double *)malloc(15 * numCells * sizeof(double));
  y_data  = (double *)malloc(15 * numCells * sizeof(double));
  xr_data = (double *)malloc(15 * numCells * sizeof(double));
  yr_data = (double *)malloc(15 * numCells * sizeof(double));
  xs_data = (double *)malloc(15 * numCells * sizeof(double));
  ys_data = (double *)malloc(15 * numCells * sizeof(double));
  rx_data = (double *)malloc(15 * numCells * sizeof(double));
  ry_data = (double *)malloc(15 * numCells * sizeof(double));
  sx_data = (double *)malloc(15 * numCells * sizeof(double));
  sy_data = (double *)malloc(15 * numCells * sizeof(double));
  nx_data = (double *)malloc(3 * 5 * numCells * sizeof(double));
  ny_data = (double *)malloc(3 * 5 * numCells * sizeof(double));
  fscale_data = (double *)malloc(3 * 5 * numCells * sizeof(double));
  for(int i = 0; i < 3; i++) {
    Q_data[i] = (double *)malloc(15 * numCells * sizeof(double));
  }
  for(int i = 0; i < 4; i++) {
    F_data[i] = (double *)malloc(15 * numCells * sizeof(double));
  }

  // Initialise OP2
  // Declare OP2 sets
  nodes  = op_decl_set(numNodes, "nodes");
  cells  = op_decl_set(numCells, "cells");
  edges  = op_decl_set(numEdges, "edges");
  bedges = op_decl_set(numBoundaryEdges, "bedges");

  // Declare OP2 maps
  cell2nodes  = op_decl_map(cells, nodes, 3, cgnsCells, "cell2nodes");
  edge2nodes  = op_decl_map(edges, nodes, 2, edge2node_data, "edge2nodes");
  edge2cells  = op_decl_map(edges, cells, 2, edge2cell_data, "edge2cells");
  bedge2nodes = op_decl_map(bedges, nodes, 2, bedge2node_data, "bedge2nodes");
  bedge2cells = op_decl_map(bedges, cells, 1, bedge2cell_data, "bedge2cells");

  // Declare OP2 datasets
    // Structure: {x, y}
  node_coords = op_decl_dat(nodes, 2, "double", coords, "node_coords");
    // Coords of nodes per cell
  nodeX = op_decl_dat(cells, 3, "double", nodeX_data, "nodeX");
  nodeY = op_decl_dat(cells, 3, "double", nodeY_data, "nodeY");
    // The x and y coordinates of all the solution points in a cell
  x = op_decl_dat(cells, 15, "double", x_data, "x");
  y = op_decl_dat(cells, 15, "double", y_data, "y");
    // Geometric factors that relate to mapping between global and local (cell) coordinates
  xr = op_decl_dat(cells, 15, "double", xr_data, "xr");
  yr = op_decl_dat(cells, 15, "double", yr_data, "yr");
  xs = op_decl_dat(cells, 15, "double", xs_data, "xs");
  ys = op_decl_dat(cells, 15, "double", ys_data, "ys");
  rx = op_decl_dat(cells, 15, "double", rx_data, "rx");
  ry = op_decl_dat(cells, 15, "double", ry_data, "ry");
  sx = op_decl_dat(cells, 15, "double", sx_data, "sx");
  sy = op_decl_dat(cells, 15, "double", sy_data, "sy");
    // Normals for each cell (calculated for each node on each edge, nodes can appear on multiple edges)
  nx = op_decl_dat(cells, 3 * 5, "double", nx_data, "nx");
  ny = op_decl_dat(cells, 3 * 5, "double", ny_data, "ny");
    // surface Jacobian / Jacobian (used when lifting the boundary fluxes)
  fscale     = op_decl_dat(cells, 3 * 5, "double", fscale_data, "fscale");
  bedge_type = op_decl_dat(bedges, 1, "int", bedge_type_data, "bedge_type");
  edgeNum    = op_decl_dat(edges, 2, "int", edgeNum_data, "edgeNum");
  bedgeNum   = op_decl_dat(bedges, 1, "int", bedgeNum_data, "bedgeNum");
  for(int i = 0; i < 3; i++) {
    string Qname = "Q" + to_string(i);
    Q[i] = op_decl_dat(cells, 15, "double", Q_data[i], Qname.c_str());
  }
  for(int i = 0; i < 4; i++) {
    string Fname = "F" + to_string(i);
    F[i] = op_decl_dat(cells, 15, "double", F_data[i], Fname.c_str());
  }
}
