#ifndef __INS_DATA_H
#define __INS_DATA_H

#include "op_seq.h"

#include "dg_global_constants.h"
#include "dg_mesh.h"

class INSData {
public:
  INSData(DGMesh *m);
  void init();

  // OP2 stuff
  op_dat Q[2][2], F[4], N[2][2], flux[2], QT[2], QTT[2], gQ[2], new_order;
  op_dat divVelT, curlVel, gradCurlVel[2], dPdN[2], pRHS, p, dpdx, dpdy, gP;
  op_dat pFluxX, pFluxY, visRHS[2], prBC, visBC[2], vorticity, save_temp;
  op_dat rho, mu, visTmp[2], gN[2], gGradCurl[2], gRho;
  op_dat tmp_dg_np[10], tmp_dg_g_np[5];

  DGMesh *mesh;
};

#endif
