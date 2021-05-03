#ifndef __INS_DATA_H
#define __INS_DATA_H

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
extern double cubV_g[46 * 15];
extern double cubW_g[46];
extern double cubVDr_g[46 * 15];
extern double cubVDs_g[46 * 15];
extern double gaussW_g[7];
extern double gFInterp0_g[7 * 15];
extern double gFInterp1_g[7 * 15];
extern double gFInterp2_g[7 * 15];
extern double gF0Dr_g[7 * 15];
extern double gF0Ds_g[7 * 15];
extern double gF1Dr_g[7 * 15];
extern double gF1Ds_g[7 * 15];
extern double gF2Dr_g[7 * 15];
extern double gF2Ds_g[7 * 15];
extern double gFInterp0R_g[7 * 15];
extern double gFInterp1R_g[7 * 15];
extern double gFInterp2R_g[7 * 15];
extern double gF0DrR_g[7 * 15];
extern double gF0DsR_g[7 * 15];
extern double gF1DrR_g[7 * 15];
extern double gF1DsR_g[7 * 15];
extern double gF2DrR_g[7 * 15];
extern double gF2DsR_g[7 * 15];
extern double lift_drag_vec[5];

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
  op_dat node_coords, nodeX, nodeY, x, y, rx, ry, sx, sy, nx,
         ny, J, sJ, fscale, bedge_type, edgeNum, bedgeNum;
  op_dat Q[2][2], exQ[2], F[4], N[2][2], flux[2], QT[2], QTT[2];
  op_dat div[4];
  op_dat divVelT, curlVel, gradCurlVel[2], dPdN[2], pRHS, pRHSex, p, dpdx, dpdy;
  op_dat visRHS[2];
  op_dat zeroBC, visBC[2];
  op_dat dQdx[2], dQdy[2];
  op_dat vorticity;
private:
  // Pointers to private memory
  double *nodeX_data;
  double *nodeY_data;
  double *x_data;
  double *y_data;
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
  double *zeroBC_data;
  double *visBC_data[2];
  double *dQdx_data[2];
  double *dQdy_data[2];
  double *vorticity_data;
};

class CubatureData {
public:
  CubatureData(INSData *dat);
  ~CubatureData();

  // mm and OP are stored in column major format
  // OP is the local stiffness matrix used by the Poisson solver
  op_dat rx, sx, ry, sy, J, mm, Dx, Dy, OP;
  op_dat temp, temp2;
  op_dat op_temps[4];

private:
  INSData *data;

  double *rx_data;
  double *sx_data;
  double *ry_data;
  double *sy_data;
  double *J_data;
  double *mm_data;
  double *Dx_data;
  double *Dy_data;
  double *OP_data;
  double *temp_data;
  double *temp2_data;
  double *op_temps_data[4];
};

class GaussData {
public:
  GaussData(INSData *dat);
  ~GaussData();

  op_dat x, y;
  op_dat rx, sx, ry, sy, sJ, nx, ny, tau, reverse;
  op_dat mDx[3], mDy[3], pDx[3], pDy[3], mD[3], pD[3];
  // OP is in column major format
  op_dat OP[3], OPf[3];
private:
  INSData *data;

  double *x_data;
  double *y_data;
  double *rx_data;
  double *sx_data;
  double *ry_data;
  double *sy_data;
  double *sJ_data;
  double *nx_data;
  double *ny_data;
  double *tau_data;
  double *mDx_data[3];
  double *mDy_data[3];
  double *pDx_data[3];
  double *pDy_data[3];
  double *mD_data[3];
  double *pD_data[3];
  double *OP_data[3];
  double *OPf_data[3];

  int *reverse_data;
};

#endif
