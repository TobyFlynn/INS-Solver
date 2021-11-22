#ifndef __INS_DATA_H
#define __INS_DATA_H

#include "op_seq.h"

#include "dg_global_constants.h"
#include "dg_mesh.h"

class INSData {
public:
  INSData(DGMesh *m);
  ~INSData();
  void init();

  // OP2 stuff
  op_dat Q[2][2], F[4], N[2][2], flux[2], QT[2], QTT[2], gQ[2];
  op_dat divVelT, curlVel, gradCurlVel[2], dPdN[2], pRHS, p, dpdx, dpdy, gP;
  op_dat pFluxX, pFluxY, visRHS[2], prBC, visBC[2], vorticity, save_temp;
  op_dat tmp_dg_np[4], tmp_dg_g_np[4];

  DGMesh *mesh;
private:
  // Pointers to private memory
  double *Q_data[2][2], *N_data[2][2], *QT_data[2], *QTT_data[2], *dPdN_data[2];
  double *p_data, *prBC_data, *vorticity_data, *save_temp_data;
  double *tmp_dg_np_data[4], *tmp_dg_g_np_data[4];
};

class CubatureData {
public:
  CubatureData(DGMesh *m, INSData *dat);
  ~CubatureData();
  void init();

  // mm and OP are stored in column major format
  // OP is the local stiffness matrix used by the Poisson solver
  op_dat Dx, Dy, OP;
  op_dat temp, temp2;

private:
  INSData *data;
  DGMesh *mesh;

  double *Dx_data;
  double *Dy_data;
  double *OP_data;
  double *temp_data;
  double *temp2_data;
};

class GaussData {
public:
  GaussData(DGMesh *m, INSData *dat);
  ~GaussData();
  void init();

  op_dat rx, sx, ry, sy, tau, reverse;
  op_dat mDx[3], mDy[3], pDx[3], pDy[3], mD[3], pD[3];
  // OP is in column major format
  op_dat OP[3], OPf[3];
private:
  INSData *data;
  DGMesh *mesh;

  double *rx_data;
  double *sx_data;
  double *ry_data;
  double *sy_data;
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
