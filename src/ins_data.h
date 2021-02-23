#ifndef __NS_DATA_H
#define __NS_DATA_H

#include "op_seq.h"

extern double gam;
extern double mu;
extern double nu;
extern double bc_mach;
extern double bc_alpha;
extern double bc_p;
extern double bc_u;
extern double bc_v;
extern int FMASK[15];
extern double ic_u;
extern double ic_v;

class INSData {
public:
  INSData();
  ~INSData();
  void initOP2();
  // Pointers used when loading data
  double *coords;
  int *cgnsCells;
  int *edge2node_data;
  int *edge2cell_data;
  int *bedge2node_data;
  int *bedge2cell_data;
  int *bedge_type_data;
  int *edgeNum_data;
  int *bedgeNum_data;
  int numNodes, numCells, numEdges, numBoundaryEdges;
  // OP2 stuff
  op_set nodes, cells, edges, bedges;
  op_map cell2nodes, edge2nodes, edge2cells, bedge2nodes, bedge2cells;
  op_dat node_coords, nodeX, nodeY, x, y, xr, yr, xs, ys, rx, ry, sx, sy, nx,
         ny, J, sJ, fscale, bedge_type, edgeNum, bedgeNum;
  op_dat Q[2][2], exQ[2], F[4], N[2][2], flux[2], QT[2], QTT[2];
  op_dat div[4];
  op_dat divVelT, curlVel, gradCurlVel[2], dPdN[2], pRHS, pRHSex, p, dpdx, dpdy;
  op_dat visRHS[2];
  op_dat dirichletBC, neumannBCx, neumannBCy;
private:
  // Pointers to private memory
  double *nodeX_data;
  double *nodeY_data;
  double *x_data;
  double *y_data;
  double *xr_data;
  double *yr_data;
  double *xs_data;
  double *ys_data;
  double *rx_data;
  double *ry_data;
  double *sx_data;
  double *sy_data;
  double *nx_data;
  double *ny_data;
  double *J_data;
  double *sJ_data;
  double *fscale_data;
  double *Q_data[2][2];
  double *exQ_data[2];
  double *F_data[4];
  double *N_data[2][2];
  double *flux_data[2];
  double *QT_data[2];
  double *QTT_data[2];
  double *div_data[4];
  double *divVelT_data;
  double *curlVel_data;
  double *gradCurlVel_data[2];
  double *dPdN_data[2];
  double *pRHS_data;
  double *pRHSex_data;
  double *p_data;
  double *dpdx_data;
  double *dpdy_data;
  double *visRHS_data[2];
  double *dirichletBC_data;
  double *neumannBCx_data;
  double *neumannBCy_data;
};

#endif
