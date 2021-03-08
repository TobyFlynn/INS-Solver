#include "ins_data.h"

#include "op_seq.h"

#include <string>
#include <memory>

#include "blas_calls.h"

#include "kernels/init_cubature_grad.h"
#include "kernels/init_cubature.h"
#include "kernels/init_cubature_OP.h"
#include "kernels/init_gauss_grad.h"
#include "kernels/gauss_grad_faces.h"
#include "kernels/init_gauss.h"

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
  op_decl_const(46, "double", cubW);
  op_decl_const(46*15, "double", cubV);
  op_decl_const(46*15, "double", cubVDr);
  op_decl_const(46*15, "double", cubVDs);
  op_decl_const(7*15, "double", gF0Dr);
  op_decl_const(7*15, "double", gF0Ds);
  op_decl_const(7*15, "double", gF1Dr);
  op_decl_const(7*15, "double", gF1Ds);
  op_decl_const(7*15, "double", gF2Dr);
  op_decl_const(7*15, "double", gF2Ds);
}

CubatureData::CubatureData(INSData *dat) {
  data = dat;

  rx_data    = (double *)malloc(46 * data->numCells * sizeof(double));
  sx_data    = (double *)malloc(46 * data->numCells * sizeof(double));
  ry_data    = (double *)malloc(46 * data->numCells * sizeof(double));
  sy_data    = (double *)malloc(46 * data->numCells * sizeof(double));
  J_data     = (double *)malloc(46 * data->numCells * sizeof(double));
  mm_data    = (double *)malloc(15 * 15 * data->numCells * sizeof(double));
  Dx_data    = (double *)malloc(46 * 15 * data->numCells * sizeof(double));
  Dy_data    = (double *)malloc(46 * 15 * data->numCells * sizeof(double));
  OP_data    = (double *)malloc(15 * 15 * data->numCells * sizeof(double));
  temp_data  = (double *)malloc(46 * 15 * data->numCells * sizeof(double));
  temp2_data = (double *)malloc(46 * 15 * data->numCells * sizeof(double));

  rx    = op_decl_dat(data->cells, 46, "double", rx_data, "cub-rx");
  sx    = op_decl_dat(data->cells, 46, "double", sx_data, "cub-sx");
  ry    = op_decl_dat(data->cells, 46, "double", ry_data, "cub-ry");
  sy    = op_decl_dat(data->cells, 46, "double", sy_data, "cub-sy");
  J     = op_decl_dat(data->cells, 46, "double", J_data, "cub-J");
  mm    = op_decl_dat(data->cells, 15 * 15, "double", mm_data, "cub-mm");
  Dx    = op_decl_dat(data->cells, 46 * 15, "double", Dx_data, "cub-Dx");
  Dy    = op_decl_dat(data->cells, 46 * 15, "double", Dy_data, "cub-Dy");
  OP    = op_decl_dat(data->cells, 15 * 15, "double", OP_data, "cub-OP");
  temp  = op_decl_dat(data->cells, 46 * 15, "double", temp_data, "cub-temp");
  temp2 = op_decl_dat(data->cells, 46 * 15, "double", temp2_data, "cub-temp2");

  // Initialise these values
  init_cubature_grad_blas(data, this);

  op_par_loop(init_cubature_grad, "init_cubature_grad", data->cells,
              op_arg_dat(rx, -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(sx, -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(ry, -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(sy, -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(Dx, -1, OP_ID, 46 * 15, "double", OP_WRITE),
              op_arg_dat(Dy, -1, OP_ID, 46 * 15, "double", OP_WRITE));

  init_cubature_blas(data, this);

  op_par_loop(init_cubature, "init_cubature", data->cells,
              op_arg_dat(rx,   -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(sx,   -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(ry,   -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(sy,   -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(J,    -1, OP_ID, 46, "double", OP_WRITE),
              op_arg_dat(temp, -1, OP_ID, 46 * 15, "double", OP_WRITE));

  cubature_mm_blas(data, this);

  op_par_loop(init_cubature_OP, "init_cubature_OP", data->cells,
              op_arg_dat(J,     -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(Dx,    -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(Dy,    -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(temp,  -1, OP_ID, 46 * 15, "double", OP_WRITE),
              op_arg_dat(temp2, -1, OP_ID, 46 * 15, "double", OP_WRITE));

  cubature_op_blas(data, this);
}

CubatureData::~CubatureData() {
  free(rx_data);
  free(sx_data);
  free(ry_data);
  free(sy_data);
  free(J_data);
  free(mm_data);
  free(Dx_data);
  free(Dy_data);
  free(OP_data);
  free(temp_data);
}

GaussData::GaussData(INSData *dat) {
  data = dat;

  rx_data = (double *)malloc(21 * data->numCells * sizeof(double));
  sx_data = (double *)malloc(21 * data->numCells * sizeof(double));
  ry_data = (double *)malloc(21 * data->numCells * sizeof(double));
  sy_data = (double *)malloc(21 * data->numCells * sizeof(double));
  sJ_data = (double *)malloc(21 * data->numCells * sizeof(double));
  nx_data = (double *)malloc(21 * data->numCells * sizeof(double));
  ny_data = (double *)malloc(21 * data->numCells * sizeof(double));
  for(int i = 0; i < 3; i++) {
    mDx_data[i] = (double *)malloc(7 * 15 * data->numCells * sizeof(double));
    mDy_data[i] = (double *)malloc(7 * 15 * data->numCells * sizeof(double));
    pDx_data[i] = (double *)calloc(7 * 15 * data->numCells, sizeof(double));
    pDy_data[i] = (double *)calloc(7 * 15 * data->numCells, sizeof(double));
  }

  rx = op_decl_dat(data->cells, 21, "double", rx_data, "gauss-rx");
  sx = op_decl_dat(data->cells, 21, "double", sx_data, "gauss-sx");
  ry = op_decl_dat(data->cells, 21, "double", ry_data, "gauss-ry");
  sy = op_decl_dat(data->cells, 21, "double", sy_data, "gauss-sy");
  sJ = op_decl_dat(data->cells, 21, "double", sJ_data, "gauss-sJ");
  nx = op_decl_dat(data->cells, 21, "double", nx_data, "gauss-nx");
  ny = op_decl_dat(data->cells, 21, "double", ny_data, "gauss-ny");
  for(int i = 0; i < 3; i++) {
    string name = "mDx" + to_string(i);
    mDx[i] = op_decl_dat(data->cells, 7 * 15, "double", mDx_data[i], name.c_str());
    name = "mDy" + to_string(i);
    mDy[i] = op_decl_dat(data->cells, 7 * 15, "double", mDy_data[i], name.c_str());
    name = "pDx" + to_string(i);
    pDx[i] = op_decl_dat(data->cells, 7 * 15, "double", pDx_data[i], name.c_str());
    name = "pDy" + to_string(i);
    pDy[i] = op_decl_dat(data->cells, 7 * 15, "double", pDy_data[i], name.c_str());
  }

  init_gauss_grad_blas(data, this);

  op_par_loop(init_gauss_grad, "init_gauss_grad", data->cells,
              op_arg_dat(rx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(ry, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sy, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(mDx[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[2], -1, OP_ID, 7 * 15, "double", OP_WRITE));

  op_par_loop(gauss_grad_faces, "gauss_grad_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(mDx[0], -2, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[0], -2, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[1], -2, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[1], -2, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[2], -2, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[2], -2, data->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(pDx[0], -2, data->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[0], -2, data->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDx[1], -2, data->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[1], -2, data->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDx[2], -2, data->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[2], -2, data->edge2cells, 7 * 15, "double", OP_INC));

  init_gauss_blas(data, this);

  op_par_loop(init_gauss, "init_gauss", data->cells,
              op_arg_dat(rx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(ry, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sy, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(nx, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(ny, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(sJ, -1, OP_ID, 21, "double", OP_WRITE));
}

GaussData::~GaussData() {
  free(rx_data);
  free(sx_data);
  free(ry_data);
  free(sy_data);
  free(sJ_data);
  free(nx_data);
  free(ny_data);
  for(int i = 0; i < 3; i++) {
    free(mDx_data[i]);
    free(mDy_data[i]);
    free(pDx_data[i]);
    free(pDy_data[i]);
  }
}
