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
  op_dat Q[2][2], F[4], N[2][2], flux[2], QT[2], QTT[2], gQ[4], new_order;
  op_dat divVelT, curlVel, gradCurlVel[2], dPdN[2], pRHS, p, dpdx, dpdy, gP;
  op_dat pFluxX, pFluxY, visRHS[2], prBC, visBC[2], vorticity, save_temp;
  op_dat rho, mu;
  op_dat tmp_dg_np[10], tmp_dg_g_np[6];

  DGMesh *mesh;
private:
  // Pointers to private memory
  double *Q_data[2][2], *N_data[2][2], *QT_data[2], *QTT_data[2], *dPdN_data[2];
  double *p_data, *prBC_data, *vorticity_data, *save_temp_data;
  double *rho_data, *mu_data;
  double *tmp_dg_np_data[10], *tmp_dg_g_np_data[6];
  int *new_order_data;
};

#endif
