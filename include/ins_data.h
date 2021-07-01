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

  op_dat Q[2][2], exQ[2], F[4], N[2][2], flux[2], QT[2], QTT[2];
  op_dat divVelT, curlVel, gradCurlVel[2], dPdN[2], pRHS, pRHSex, p, dpdx, dpdy;
  op_dat visRHS[2], prBC, visBC[2], dQdx[2], dQdy[2];
  op_dat vorticity, save_temp, nu, gNu, rho, pFluxX, pFluxY;

  // Cubature stuff
  op_dat Dx, Dy, cOP, temp, temp2;

  // Gauss stuff
  op_dat grx, gsx, gry, gsy, tau, reverse;
  op_dat mDx[3], mDy[3], pDx[3], pDy[3], mD[3], pD[3];
  // OP is in column major format
  op_dat gOP[3], gOPf[3];
private:
  DGMesh *mesh;

  // Pointers to private memory
  double *Q_data[2][2], *exQ_data[2], *F_data[4], *N_data[2][2];
  double *flux_data[2], *QT_data[2], *QTT_data[2];
  double *divVelT_data, *curlVel_data, *gradCurlVel_data[2];
  double *dPdN_data[2], *pRHS_data, *pRHSex_data, *p_data;
  double *dpdx_data, *dpdy_data, *prBC_data;
  double *visRHS_data[2], *visBC_data[2];
  double *dQdx_data[2], *dQdy_data[2];
  double *vorticity_data, *save_temp_data;
  double *nu_data, *gNu_data, *rho_data, *pFluxX_data, *pFluxY_data;

  // Cubature stuff
  double *Dx_data, *Dy_data, *cOP_data, *temp_data, *temp2_data;

  // Gauss stuff
  double *grx_data, *gsx_data, *gry_data, *gsy_data, *tau_data;
  int *reverse_data;
  double *mDx_data[3], *mDy_data[3];
  double *pDx_data[3], *pDy_data[3];
  double *mD_data[3], *pD_data[3];
  double *gOP_data[3], *gOPf_data[3];
};

#endif
