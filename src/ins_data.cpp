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
  free(J_data);
  free(sJ_data);
  free(fscale_data);
  for(int i = 0; i < 4; i++) {
    free(F_data[i]);
    free(div_data[i]);
  }
  for(int i = 0; i < 2; i++) {
    free(Q_data[0][i]);
    free(Q_data[1][i]);
    free(QT_data[i]);
    free(QTT_data[i]);

    free(N_data[0][i]);
    free(N_data[1][i]);
    free(exQ_data[i]);
    free(flux_data[i]);
    free(gradCurlVel_data[i]);
    free(dPdN_data[i]);
    free(visRHS_data[i]);
    free(neumannBCx_data[i]);
    free(neumannBCy_data[i]);
  }
  free(divVelT_data);
  free(curlVel_data);
  free(pRHS_data);
  free(pRHSex_data);
  free(p_data);
  free(dpdx_data);
  free(dpdy_data);
  free(dirichletBC_data);
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
  nx_data = (double *)malloc(15 * numCells * sizeof(double));
  ny_data = (double *)malloc(15 * numCells * sizeof(double));
  J_data  = (double *)malloc(15 * numCells * sizeof(double));
  sJ_data = (double *)malloc(15 * numCells * sizeof(double));
  fscale_data = (double *)malloc(15 * numCells * sizeof(double));
  for(int i = 0; i < 4; i++) {
    F_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    div_data[i] = (double *)malloc(15 * numCells * sizeof(double));
  }
  for(int i = 0; i < 2; i++) {
    Q_data[0][i] = (double *)malloc(15 * numCells * sizeof(double));
    Q_data[1][i] = (double *)malloc(15 * numCells * sizeof(double));
    QT_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    QTT_data[i] = (double *)malloc(15 * numCells * sizeof(double));

    N_data[0][i] = (double *)malloc(15 * numCells * sizeof(double));
    N_data[1][i] = (double *)malloc(15 * numCells * sizeof(double));
    exQ_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    flux_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    gradCurlVel_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    dPdN_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    visRHS_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    neumannBCx_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    neumannBCy_data[i] = (double *)malloc(15 * numCells * sizeof(double));
  }
  divVelT_data = (double *)malloc(15 * numCells * sizeof(double));
  curlVel_data = (double *)malloc(15 * numCells * sizeof(double));
  pRHS_data    = (double *)malloc(15 * numCells * sizeof(double));
  pRHSex_data  = (double *)malloc(15 * numCells * sizeof(double));
  p_data       = (double *)malloc(15 * numCells * sizeof(double));
  dpdx_data    = (double *)malloc(15 * numCells * sizeof(double));
  dpdy_data    = (double *)malloc(15 * numCells * sizeof(double));
  dirichletBC_data = (double *)malloc(15 * numCells * sizeof(double));

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
  nx = op_decl_dat(cells, 15, "double", nx_data, "nx");
  ny = op_decl_dat(cells, 15, "double", ny_data, "ny");
    // surface Jacobian / Jacobian (used when lifting the boundary fluxes)
  J          = op_decl_dat(cells, 15, "double", J_data, "J");
  sJ         = op_decl_dat(cells, 15, "double", sJ_data, "sJ");
  fscale     = op_decl_dat(cells, 15, "double", fscale_data, "fscale");
  bedge_type = op_decl_dat(bedges, 1, "int", bedge_type_data, "bedge_type");
  edgeNum    = op_decl_dat(edges, 2, "int", edgeNum_data, "edgeNum");
  bedgeNum   = op_decl_dat(bedges, 1, "int", bedgeNum_data, "bedgeNum");
  for(int i = 0; i < 4; i++) {
    string Fname = "F" + to_string(i);
    F[i] = op_decl_dat(cells, 15, "double", F_data[i], Fname.c_str());
    string divname = "div" + to_string(i);
    div[i] = op_decl_dat(cells, 15, "double", div_data[i], divname.c_str());
  }
  for(int i = 0; i < 2; i++) {
    string Qname = "Q0" + to_string(i);
    Q[0][i] = op_decl_dat(cells, 15, "double", Q_data[0][i], Qname.c_str());
    Qname = "Q1" + to_string(i);
    Q[1][i] = op_decl_dat(cells, 15, "double", Q_data[1][i], Qname.c_str());
    Qname = "QT" + to_string(i);
    QT[i] = op_decl_dat(cells, 15, "double", QT_data[i], Qname.c_str());
    Qname = "QTT" + to_string(i);
    QTT[i] = op_decl_dat(cells, 15, "double", QTT_data[i], Qname.c_str());

    string Nname = "N0" + to_string(i);
    N[0][i] = op_decl_dat(cells, 15, "double", N_data[0][i], Nname.c_str());
    Nname = "N1" + to_string(i);
    N[1][i] = op_decl_dat(cells, 15, "double", N_data[1][i], Nname.c_str());
    string exQname = "exQ" + to_string(i);
    exQ[i] = op_decl_dat(cells, 15, "double", exQ_data[i], exQname.c_str());
    string fluxname = "flux" + to_string(i);
    flux[i] = op_decl_dat(cells, 15, "double", flux_data[i], fluxname.c_str());
    string gradCurlVelname = "gradCurlVel" + to_string(i);
    gradCurlVel[i] = op_decl_dat(cells, 15, "double", gradCurlVel_data[i], gradCurlVelname.c_str());
    string dPdNname = "dPdN" + to_string(i);
    dPdN[i] = op_decl_dat(cells, 15, "double", dPdN_data[i], dPdNname.c_str());
    string visRHSname = "visRHS" + to_string(i);
    visRHS[i] = op_decl_dat(cells, 15, "double", visRHS_data[i], visRHSname.c_str());
    neumannBCx[i] = op_decl_dat(cells, 15, "double", neumannBCx_data[i], "neumannBCx");
    neumannBCy[i] = op_decl_dat(cells, 15, "double", neumannBCy_data[i], "neumannBCy");
  }
  divVelT = op_decl_dat(cells, 15, "double", divVelT_data, "divVelT");
  curlVel = op_decl_dat(cells, 15, "double", curlVel_data, "curlVel");
  pRHS    = op_decl_dat(cells, 15, "double", pRHS_data, "pRHS");
  pRHSex  = op_decl_dat(cells, 15, "double", pRHSex_data, "pRHSex");
  p       = op_decl_dat(cells, 15, "double", p_data, "p");
  dpdx    = op_decl_dat(cells, 15, "double", dpdx_data, "dpdx");
  dpdy    = op_decl_dat(cells, 15, "double", dpdy_data, "dpdy");
  dirichletBC = op_decl_dat(cells, 15, "double", dirichletBC_data, "dirichletBC");

  op_decl_const(1, "double", &gam);
  op_decl_const(1, "double", &mu);
  op_decl_const(1, "double", &nu);
  op_decl_const(1, "double", &bc_mach);
  op_decl_const(1, "double", &bc_alpha);
  op_decl_const(1, "double", &bc_p);
  op_decl_const(1, "double", &bc_u);
  op_decl_const(1, "double", &bc_v);
  op_decl_const(15, "int", FMASK);
  op_decl_const(1, "double", &ic_u);
  op_decl_const(1, "double", &ic_v);
}
