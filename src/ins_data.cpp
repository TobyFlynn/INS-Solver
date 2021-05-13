#include "ins_data.h"

#include "op_seq.h"

#include <string>
#include <memory>

#include "blas_calls.h"
#include "load_mesh.h"

#include "kernels/init_grid.h"
#include "kernels/init_nodes.h"

#include "kernels/init_cubature_grad.h"
#include "kernels/init_cubature.h"
#include "kernels/init_cubature_OP.h"

#include "kernels/gauss_reverse.h"
#include "kernels/gauss_tau.h"
#include "kernels/gauss_tau_bc.h"
#include "kernels/init_gauss_grad.h"
#include "kernels/gauss_grad_faces.h"
#include "kernels/init_gauss_grad2.h"
#include "kernels/init_gauss.h"
#include "kernels/gauss_op.h"
#include "kernels/gauss_gfi_faces.h"
#include "kernels/init_gauss_grad_neighbour.h"

using namespace std;

INSData::INSData(std::string filename) {
  #ifdef POISSON_TEST
  // Lamda used to identify the type of boundary edges
  auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
    if(y1 == y2 && y1 > 0.5) {
      // Neumann BC y = 1
      return 2;
    } else if(y1 == y2 && y1 < 0.5) {
      // Neumann BC y = 0
      return 3;
    } else if(x1 < 0.5){
      // Dirichlet BC x = 0
      return 0;
    } else {
      // Dirichlet BC x = 1
      return 1;
    }
  };
  // auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
  //   if(y1 == y2 && y1 > 0.5) {
  //     // Neumann BC y = 1
  //     return 1;
  //   } else if(y1 == y2 && y1 < 0.5) {
  //     // Neumann BC y = 0
  //     return 1;
  //   } else if(x1 < 0.5){
  //     // Dirichlet BC x = 0
  //     return 1;
  //   } else {
  //     // Dirichlet BC x = 1
  //     return 0;
  //   }
  // };
  #else
  // Lamda used to identify the type of boundary edges
  auto bcNum = [](double x1, double x2, double y1, double y2) -> int {
    if(x1 == 0.0 && x2 == 0.0) {
      // Inflow
      return 0;
    } else if(x1 == 2.2 && x2 == 2.2) {
      // Outflow
      return 1;
    } else if(x1 > 0.1 && x2 > 0.1 && x1 < 0.3 && x2 < 0.3
              && y1 > 0.1 && y2 > 0.1 && y1 < 0.3 && y2 < 0.3) {
      // Cylinder Wall
      return 2;
    } else {
      // Top/Bottom Wall
      return 3;
    }
  };
  #endif

  load_mesh(filename.c_str(), this, bcNum);

  // Initialise memory
  nodeX_data  = (double*)calloc(3 * numCells, sizeof(double));
  nodeY_data  = (double*)calloc(3 * numCells, sizeof(double));
  x_data      = (double *)calloc(15 * numCells, sizeof(double));
  y_data      = (double *)calloc(15 * numCells, sizeof(double));
  rx_data     = (double *)calloc(15 * numCells, sizeof(double));
  ry_data     = (double *)calloc(15 * numCells, sizeof(double));
  sx_data     = (double *)calloc(15 * numCells, sizeof(double));
  sy_data     = (double *)calloc(15 * numCells, sizeof(double));
  nx_data     = (double *)calloc(15 * numCells, sizeof(double));
  ny_data     = (double *)calloc(15 * numCells, sizeof(double));
  J_data      = (double *)calloc(15 * numCells, sizeof(double));
  sJ_data     = (double *)calloc(15 * numCells, sizeof(double));
  fscale_data = (double *)calloc(15 * numCells, sizeof(double));
  for(int i = 0; i < 4; i++) {
    F_data[i]   = (double *)calloc(15 * numCells, sizeof(double));
    div_data[i] = (double *)calloc(15 * numCells, sizeof(double));
  }
  for(int i = 0; i < 2; i++) {
    Q_data[0][i]   = (double *)calloc(15 * numCells, sizeof(double));
    Q_data[1][i]   = (double *)calloc(15 * numCells, sizeof(double));
    QT_data[i]     = (double *)calloc(15 * numCells, sizeof(double));
    QTT_data[i]    = (double *)calloc(15 * numCells, sizeof(double));
    N_data[0][i]   = (double *)calloc(15 * numCells, sizeof(double));
    N_data[1][i]   = (double *)calloc(15 * numCells, sizeof(double));
    exQ_data[i]    = (double *)calloc(15 * numCells, sizeof(double));
    flux_data[i]   = (double *)calloc(15 * numCells, sizeof(double));
    dPdN_data[i]   = (double *)calloc(15 * numCells, sizeof(double));
    visRHS_data[i] = (double *)calloc(15 * numCells, sizeof(double));
    visBC_data[i]  = (double *)calloc(21 * numCells, sizeof(double));
    dQdx_data[i]   = (double *)calloc(15 * numCells, sizeof(double));
    dQdy_data[i]   = (double *)calloc(15 * numCells, sizeof(double));
    gradCurlVel_data[i] = (double *)calloc(15 * numCells, sizeof(double));
  }
  divVelT_data   = (double *)calloc(15 * numCells, sizeof(double));
  curlVel_data   = (double *)calloc(15 * numCells, sizeof(double));
  pRHS_data      = (double *)calloc(15 * numCells, sizeof(double));
  pRHSex_data    = (double *)calloc(15 * numCells, sizeof(double));
  p_data         = (double *)calloc(15 * numCells, sizeof(double));
  dpdx_data      = (double *)calloc(15 * numCells, sizeof(double));
  dpdy_data      = (double *)calloc(15 * numCells, sizeof(double));
  zeroBC_data    = (double *)calloc(21 * numCells, sizeof(double));
  vorticity_data = (double *)calloc(15 * numCells, sizeof(double));

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
    string dQdxname = "dQdx" + to_string(i);
    dQdx[i] = op_decl_dat(cells, 15, "double", dQdx_data[i], dQdxname.c_str());
    string dQdyname = "dQdy" + to_string(i);
    dQdy[i] = op_decl_dat(cells, 15, "double", dQdy_data[i], dQdyname.c_str());
    string visBCname = "visBC" + to_string(i);
    visBC[i] = op_decl_dat(cells, 21, "double", visBC_data[i], visBCname.c_str());
  }
  divVelT   = op_decl_dat(cells, 15, "double", divVelT_data, "divVelT");
  curlVel   = op_decl_dat(cells, 15, "double", curlVel_data, "curlVel");
  pRHS      = op_decl_dat(cells, 15, "double", pRHS_data, "pRHS");
  pRHSex    = op_decl_dat(cells, 15, "double", pRHSex_data, "pRHSex");
  p         = op_decl_dat(cells, 15, "double", p_data, "p");
  dpdx      = op_decl_dat(cells, 15, "double", dpdx_data, "dpdx");
  dpdy      = op_decl_dat(cells, 15, "double", dpdy_data, "dpdy");
  zeroBC    = op_decl_dat(cells, 21, "double", zeroBC_data, "zeroBC");
  vorticity = op_decl_dat(cells, 15, "double", vorticity_data, "vorticity");

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
  op_decl_const(46, "double", cubW_g);
  op_decl_const(46*15, "double", cubV_g);
  op_decl_const(46*15, "double", cubVDr_g);
  op_decl_const(46*15, "double", cubVDs_g);
  op_decl_const(7*15, "double", gF0Dr_g);
  op_decl_const(7*15, "double", gF0Ds_g);
  op_decl_const(7*15, "double", gF1Dr_g);
  op_decl_const(7*15, "double", gF1Ds_g);
  op_decl_const(7*15, "double", gF2Dr_g);
  op_decl_const(7*15, "double", gF2Ds_g);
  op_decl_const(7, "double", gaussW_g);
  op_decl_const(7*15, "double", gFInterp0_g);
  op_decl_const(7*15, "double", gFInterp1_g);
  op_decl_const(7*15, "double", gFInterp2_g);
  op_decl_const(7*15, "double", gF0DrR_g);
  op_decl_const(7*15, "double", gF0DsR_g);
  op_decl_const(7*15, "double", gF1DrR_g);
  op_decl_const(7*15, "double", gF1DsR_g);
  op_decl_const(7*15, "double", gF2DrR_g);
  op_decl_const(7*15, "double", gF2DsR_g);
  op_decl_const(7*15, "double", gFInterp0R_g);
  op_decl_const(7*15, "double", gFInterp1R_g);
  op_decl_const(7*15, "double", gFInterp2R_g);
  op_decl_const(5, "double", lift_drag_vec);
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
    free(visBC_data[i]);
    free(dQdx_data[i]);
    free(dQdy_data[i]);
  }
  free(divVelT_data);
  free(curlVel_data);
  free(pRHS_data);
  free(pRHSex_data);
  free(p_data);
  free(dpdx_data);
  free(dpdy_data);
  free(zeroBC_data);
  free(vorticity_data);
}

void INSData::init() {
  op_par_loop(init_nodes, "init_nodes", cells,
              op_arg_dat(node_coords, -3, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(nodeX, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(nodeY, -1, OP_ID, 3, "double", OP_WRITE));

  // Calculate geometric factors
  init_grid_blas(this);

  op_par_loop(init_grid, "init_grid", cells,
              op_arg_dat(rx, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(ry, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(sx, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(sy, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(nx, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(ny, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(J,  -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(sJ, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(fscale, -1, OP_ID, 15, "double", OP_WRITE));
}

CubatureData::CubatureData(INSData *dat) {
  data = dat;

  rx_data    = (double *)calloc(46 * data->numCells, sizeof(double));
  sx_data    = (double *)calloc(46 * data->numCells, sizeof(double));
  ry_data    = (double *)calloc(46 * data->numCells, sizeof(double));
  sy_data    = (double *)calloc(46 * data->numCells, sizeof(double));
  J_data     = (double *)calloc(46 * data->numCells, sizeof(double));
  mm_data    = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  Dx_data    = (double *)calloc(46 * 15 * data->numCells, sizeof(double));
  Dy_data    = (double *)calloc(46 * 15 * data->numCells, sizeof(double));
  OP_data    = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  temp_data  = (double *)calloc(46 * 15 * data->numCells, sizeof(double));
  temp2_data = (double *)calloc(46 * 15 * data->numCells, sizeof(double));

  op_temps_data[0] = (double *)calloc(46 * data->numCells, sizeof(double));
  op_temps_data[1] = (double *)calloc(46 * data->numCells, sizeof(double));
  op_temps_data[2] = (double *)calloc(46 * data->numCells, sizeof(double));
  op_temps_data[3] = (double *)calloc(46 * data->numCells, sizeof(double));

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

  op_temps[0] = op_decl_dat(data->cells, 46, "double", op_temps_data[0], "cub-op-temp0");
  op_temps[1] = op_decl_dat(data->cells, 46, "double", op_temps_data[1], "cub-op-temp1");
  op_temps[2] = op_decl_dat(data->cells, 46, "double", op_temps_data[2], "cub-op-temp2");
  op_temps[3] = op_decl_dat(data->cells, 46, "double", op_temps_data[3], "cub-op-temp3");
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
  free(temp2_data);
  free(op_temps_data[0]);
  free(op_temps_data[1]);
  free(op_temps_data[2]);
  free(op_temps_data[3]);
}

void CubatureData::init() {
  // Initialise geometric factors for calcuating grad matrix
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_VDR), 15, data->x, 0.0, rx);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_VDS), 15, data->x, 0.0, sx);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_VDR), 15, data->y, 0.0, ry);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_VDS), 15, data->y, 0.0, sy);

  op_par_loop(init_cubature_grad, "init_cubature_grad", data->cells,
              op_arg_dat(rx, -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(sx, -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(ry, -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(sy, -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(Dx, -1, OP_ID, 46 * 15, "double", OP_WRITE),
              op_arg_dat(Dy, -1, OP_ID, 46 * 15, "double", OP_WRITE));
  // Dx and Dy are row-major at this point

  // Calculate geometric factors for cubature volume nodes
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DR), 15, data->x, 0.0, rx);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DS), 15, data->x, 0.0, sx);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DR), 15, data->y, 0.0, ry);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DS), 15, data->y, 0.0, sy);

  op_par_loop(init_cubature, "init_cubature", data->cells,
              op_arg_dat(rx,   -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(sx,   -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(ry,   -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(sy,   -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(J,    -1, OP_ID, 46, "double", OP_WRITE),
              op_arg_dat(temp, -1, OP_ID, 46 * 15, "double", OP_WRITE));
  // Temp is in row-major at this point
  op2_gemm(false, true, 15, 15, 46, 1.0, constants->get_ptr(Constants::CUB_V), 15, temp, 15, 0.0, mm, 15);
  // mm is in col-major at this point

  // Calculate Cubature OP (contribution of Cubature points to Poisson matrix)
  op_par_loop(init_cubature_OP, "init_cubature_OP", data->cells,
              op_arg_dat(J,     -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(Dx,    -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(Dy,    -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(temp,  -1, OP_ID, 46 * 15, "double", OP_WRITE),
              op_arg_dat(temp2, -1, OP_ID, 46 * 15, "double", OP_WRITE));
  // Temp and temp2 are in row-major at this point
  op2_gemm_batch(false, true, 15, 15, 46, 1.0, Dx, 15, temp, 15, 0.0, OP, 15);
  op2_gemm_batch(false, true, 15, 15, 46, 1.0, Dy, 15, temp2, 15, 1.0, OP, 15);
  // OP is in col-major at this point
}

GaussData::GaussData(INSData *dat) {
  data = dat;

  x_data       = (double *)calloc(21 * data->numCells, sizeof(double));
  y_data       = (double *)calloc(21 * data->numCells, sizeof(double));
  rx_data      = (double *)calloc(21 * data->numCells, sizeof(double));
  sx_data      = (double *)calloc(21 * data->numCells, sizeof(double));
  ry_data      = (double *)calloc(21 * data->numCells, sizeof(double));
  sy_data      = (double *)calloc(21 * data->numCells, sizeof(double));
  sJ_data      = (double *)calloc(21 * data->numCells, sizeof(double));
  nx_data      = (double *)calloc(21 * data->numCells, sizeof(double));
  ny_data      = (double *)calloc(21 * data->numCells, sizeof(double));
  tau_data     = (double *)calloc(3 * data->numCells, sizeof(double));
  reverse_data = (int *)calloc(3 * data->numCells, sizeof(int));
  for(int i = 0; i < 3; i++) {
    mDx_data[i] = (double *)calloc(7 * 15 * data->numCells, sizeof(double));
    mDy_data[i] = (double *)calloc(7 * 15 * data->numCells, sizeof(double));
    pDx_data[i] = (double *)calloc(7 * 15 * data->numCells, sizeof(double));
    pDy_data[i] = (double *)calloc(7 * 15 * data->numCells, sizeof(double));
    mD_data[i]  = (double *)calloc(7 * 15 * data->numCells, sizeof(double));
    pD_data[i]  = (double *)calloc(7 * 15 * data->numCells, sizeof(double));
    OP_data[i]  = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
    OPf_data[i] = (double *)calloc(15 * 15 * data->numCells, sizeof(double));
  }

  x       = op_decl_dat(data->cells, 21, "double", x_data, "gauss-x");
  y       = op_decl_dat(data->cells, 21, "double", y_data, "gauss-y");
  rx      = op_decl_dat(data->cells, 21, "double", rx_data, "gauss-rx");
  sx      = op_decl_dat(data->cells, 21, "double", sx_data, "gauss-sx");
  ry      = op_decl_dat(data->cells, 21, "double", ry_data, "gauss-ry");
  sy      = op_decl_dat(data->cells, 21, "double", sy_data, "gauss-sy");
  sJ      = op_decl_dat(data->cells, 21, "double", sJ_data, "gauss-sJ");
  nx      = op_decl_dat(data->cells, 21, "double", nx_data, "gauss-nx");
  ny      = op_decl_dat(data->cells, 21, "double", ny_data, "gauss-ny");
  tau     = op_decl_dat(data->cells, 3, "double", tau_data, "gauss-tau");
  reverse = op_decl_dat(data->cells, 3, "int", reverse_data, "gauss-reverse");
  for(int i = 0; i < 3; i++) {
    string name = "mDx" + to_string(i);
    mDx[i] = op_decl_dat(data->cells, 7 * 15, "double", mDx_data[i], name.c_str());
    name = "mDy" + to_string(i);
    mDy[i] = op_decl_dat(data->cells, 7 * 15, "double", mDy_data[i], name.c_str());
    name = "pDx" + to_string(i);
    pDx[i] = op_decl_dat(data->cells, 7 * 15, "double", pDx_data[i], name.c_str());
    name = "pDy" + to_string(i);
    pDy[i] = op_decl_dat(data->cells, 7 * 15, "double", pDy_data[i], name.c_str());
    name = "mD" + to_string(i);
    mD[i] = op_decl_dat(data->cells, 7 * 15, "double", mD_data[i], name.c_str());
    name = "pD" + to_string(i);
    pD[i] = op_decl_dat(data->cells, 7 * 15, "double", pD_data[i], name.c_str());
    name = "OP" + to_string(i);
    OP[i] = op_decl_dat(data->cells, 15 * 15, "double", OP_data[i], name.c_str());
    name = "OPf" + to_string(i);
    OPf[i] = op_decl_dat(data->cells, 15 * 15, "double", OPf_data[i], name.c_str());
  }
}

GaussData::~GaussData() {
  free(x_data);
  free(y_data);
  free(rx_data);
  free(sx_data);
  free(ry_data);
  free(sy_data);
  free(sJ_data);
  free(nx_data);
  free(ny_data);
  free(tau_data);
  free(reverse_data);
  for(int i = 0; i < 3; i++) {
    free(mDx_data[i]);
    free(mDy_data[i]);
    free(pDx_data[i]);
    free(pDy_data[i]);
    free(mD_data[i]);
    free(pD_data[i]);
    free(OP_data[i]);
    free(OPf_data[i]);
  }
}

void GaussData::init() {
  // Check which edges will require matrices to be 'reverse'
  op_par_loop(gauss_reverse, "gauss_reverse", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(reverse, -2, data->edge2cells, 3, "int", OP_INC));

  // Initialise geometric factors for Gauss nodes
  init_gauss_blas(data, this);

  op_par_loop(init_gauss, "init_gauss", data->cells,
              op_arg_dat(rx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(ry, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sy, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(nx, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(ny, -1, OP_ID, 21, "double", OP_WRITE),
              op_arg_dat(sJ, -1, OP_ID, 21, "double", OP_WRITE));

  // Calculate x and y coords of Gauss nodes
  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, data->x, 0.0, x);
  op2_gemv(true, 21, 15, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, data->y, 0.0, y);

  // Calculate tau (used when constructing the Poisson matrix)
  op_par_loop(gauss_tau, "gauss_tau", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->fscale, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(tau, -2, data->edge2cells, 3, "double", OP_INC));

  op_par_loop(gauss_tau_bc, "gauss_tau_bc", data->bedges,
              op_arg_dat(data->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->fscale, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(tau, 0, data->bedge2cells, 3, "double", OP_INC));

  // Calculate geometric factors used when constructing gradient matrices
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

  // Construct gradient matrices (1 for each face)
  op_par_loop(init_gauss_grad2, "init_gauss_grad2", data->cells,
              op_arg_dat(nx, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(ny, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(mDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mD[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mD[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mD[2], -1, OP_ID, 7 * 15, "double", OP_WRITE));

  // Calculate geometric factors for grad matrices used by neighbours
  // Matrices are calculated locally, then copied to neighbour elements
  init_gauss_grad_neighbour_blas(data, this);

  op_par_loop(init_gauss_grad_neighbour, "init_gauss_grad_neighbour", data->cells,
              op_arg_dat(reverse, -1, OP_ID, 3, "int", OP_READ),
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

  // Copy x and y grad matrices to neighbours
  op_par_loop(gauss_grad_faces, "gauss_grad_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
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

  // Calculate final neighbour grad matrix
  op_par_loop(init_gauss_grad2, "init_gauss_grad2", data->cells,
              op_arg_dat(nx, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(ny, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(pDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDy[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDy[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDx[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDy[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pD[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pD[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pD[2], -1, OP_ID, 7 * 15, "double", OP_WRITE));

  // Calculate Gauss OP for each face (local contribution of face in Poisson matrix)
  // Face 0 temps: mDx, Face 1 temps: mDy, Face 2 temps: pDx
  op_par_loop(gauss_op, "gauss_op", data->cells,
              op_arg_dat(tau, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(sJ, -1, OP_ID, 21, "double", OP_READ),
              // Face 0
              op_arg_dat(mD[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              // Face 1
              op_arg_dat(mD[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              // Face 2
              op_arg_dat(mD[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDx[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pDx[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pDx[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              // Reset dats for OPf
              op_arg_dat(pDy[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pDy[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pDy[2], -1, OP_ID, 7 * 15, "double", OP_WRITE));

  op2_gemm(true, true, 15, 15, 7, 1.0, mDx[0], 7, constants->get_ptr(Constants::GAUSS_FINTERP0), 15, 0.0, OP[0], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, mDx[1], 7, mD[0], 15, 1.0, OP[0], 15);
  op2_gemm(true, true, 15, 15, 7, -1.0, mDx[2], 7, constants->get_ptr(Constants::GAUSS_FINTERP0), 15, 1.0, OP[0], 15);

  op2_gemm(true, true, 15, 15, 7, 1.0, mDy[0], 7, constants->get_ptr(Constants::GAUSS_FINTERP1), 15, 0.0, OP[1], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, mDy[1], 7, mD[1], 15, 1.0, OP[1], 15);
  op2_gemm(true, true, 15, 15, 7, -1.0, mDy[2], 7, constants->get_ptr(Constants::GAUSS_FINTERP1), 15, 1.0, OP[1], 15);

  op2_gemm(true, true, 15, 15, 7, 1.0, pDx[0], 7, constants->get_ptr(Constants::GAUSS_FINTERP2), 15, 0.0, OP[2], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, pDx[1], 7, mD[2], 15, 1.0, OP[2], 15);
  op2_gemm(true, true, 15, 15, 7, -1.0, pDx[2], 7, constants->get_ptr(Constants::GAUSS_FINTERP2), 15, 1.0, OP[2], 15);

  // Calculate Gauss OPf for each face (contribution to neighbouring element in Poisson matrix)
  op_par_loop(gauss_gfi_faces, "gauss_gfi_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(pDy[0], -2, data->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[1], -2, data->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[2], -2, data->edge2cells, 7 * 15, "double", OP_INC));

  op2_gemm_batch(true, true, 15, 15, 7, 1.0, mDx[0], 7, pDy[0], 15, 0.0, OPf[0], 15);
  op2_gemm_batch(true, true, 15, 15, 7, 1.0, mDx[1], 7, pD[0], 15, 1.0, OPf[0], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, mDx[2], 7, pDy[0], 15, 1.0, OPf[0], 15);

  op2_gemm_batch(true, true, 15, 15, 7, 1.0, mDy[0], 7, pDy[1], 15, 0.0, OPf[1], 15);
  op2_gemm_batch(true, true, 15, 15, 7, 1.0, mDy[1], 7, pD[1], 15, 1.0, OPf[1], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, mDy[2], 7, pDy[1], 15, 1.0, OPf[1], 15);

  op2_gemm_batch(true, true, 15, 15, 7, 1.0, pDx[0], 7, pDy[2], 15, 0.0, OPf[2], 15);
  op2_gemm_batch(true, true, 15, 15, 7, 1.0, pDx[1], 7, pD[2], 15, 1.0, OPf[2], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, pDx[2], 7, pDy[2], 15, 1.0, OPf[2], 15);

  // Applying the correct factors to OP and OPf is done when constructing the Poisson matrix
}
