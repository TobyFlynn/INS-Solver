#include "ins_data.h"

#include "op_seq.h"

#include <string>
#include <memory>

#include "dg_blas_calls.h"
#include "dg_compiler_defs.h"

#include "constants.h"

using namespace std;

INSData::INSData(DGMesh *m) {
  mesh = m;
  // Initialise memory
  for(int i = 0; i < 10; i++) {
    tmp_dg_np_data[i]   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  }
  for(int i = 0; i < 5; i++) {
    tmp_dg_g_np_data[i] = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  }
  for(int i = 0; i < 2; i++) {
    Q_data[0][i]   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    Q_data[1][i]   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    QT_data[i]     = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    QTT_data[i]    = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    N_data[0][i]   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    N_data[1][i]   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    dPdN_data[i]   = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  }
  p_data         = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  prBC_data      = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  vorticity_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  save_temp_data = (double *)calloc(DG_SUB_CELLS * mesh->numCells, sizeof(double));
  new_order_data = (int *)calloc(mesh->numCells, sizeof(int));
  rho_data       = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  mu_data        = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  // Declare OP2 datasets
  for(int i = 0; i < 10; i++) {
    string name    = "tmp_dg_np" + to_string(i);
    tmp_dg_np[i]   = op_decl_dat(mesh->cells, DG_NP, "double", tmp_dg_np_data[i], name.c_str());
  }
  for(int i = 0; i < 5; i++) {
    string name    = "tmp_dg_g_np" + to_string(i);
    tmp_dg_g_np[i] = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_dg_g_np_data[i], name.c_str());
  }
  for(int i = 0; i < 2; i++) {
    string name = "Q0" + to_string(i);
    Q[0][i]     = op_decl_dat(mesh->cells, DG_NP, "double", Q_data[0][i], name.c_str());
    name        = "Q1" + to_string(i);
    Q[1][i]     = op_decl_dat(mesh->cells, DG_NP, "double", Q_data[1][i], name.c_str());
    name        = "QT" + to_string(i);
    QT[i]       = op_decl_dat(mesh->cells, DG_NP, "double", QT_data[i], name.c_str());
    name        = "QTT" + to_string(i);
    QTT[i]      = op_decl_dat(mesh->cells, DG_NP, "double", QTT_data[i], name.c_str());
    name        = "N0" + to_string(i);
    N[0][i]     = op_decl_dat(mesh->cells, DG_NP, "double", N_data[0][i], name.c_str());
    name        = "N1" + to_string(i);
    N[1][i]     = op_decl_dat(mesh->cells, DG_NP, "double", N_data[1][i], name.c_str());
    name        = "dPdN" + to_string(i);
    dPdN[i]     = op_decl_dat(mesh->cells, DG_G_NP, "double", dPdN_data[i], name.c_str());
  }
  p         = op_decl_dat(mesh->cells, DG_NP, "double", p_data, "p");
  prBC      = op_decl_dat(mesh->cells, DG_G_NP, "double", prBC_data, "prBC");
  vorticity = op_decl_dat(mesh->cells, DG_NP, "double", vorticity_data, "vorticity");
  save_temp = op_decl_dat(mesh->cells, DG_SUB_CELLS, "double", save_temp_data, "save_temp");
  new_order = op_decl_dat(mesh->cells, 1, "int", new_order_data, "new_order");
  rho       = op_decl_dat(mesh->cells, DG_NP, "double", rho_data, "rho");
  mu        = op_decl_dat(mesh->cells, DG_NP, "double", mu_data, "mu");

  op_decl_const(1, "double", &reynolds);
  op_decl_const(1, "double", &mu0);
  op_decl_const(1, "double", &mu1);
  op_decl_const(1, "double", &rho0);
  op_decl_const(1, "double", &rho1);
  op_decl_const(1, "double", &nu);
  op_decl_const(1, "double", &ic_u);
  op_decl_const(1, "double", &ic_v);
  op_decl_const(3 * 5, "int", DG_CONSTANTS);
  op_decl_const(3 * 3 * DG_NPF, "int", FMASK);
  op_decl_const(3 * DG_CUB_NP, "double", cubW_g);
  op_decl_const(3 * DG_GF_NP, "double", gaussW_g);
}

INSData::~INSData() {
  for(int i = 0; i < 10; i++) {
    free(tmp_dg_np_data[i]);
  }
  for(int i = 0; i < 5; i++) {
    free(tmp_dg_g_np_data[i]);
  }
  for(int i = 0; i < 2; i++) {
    free(Q_data[0][i]);
    free(Q_data[1][i]);
    free(QT_data[i]);
    free(QTT_data[i]);
    free(N_data[0][i]);
    free(N_data[1][i]);
    free(dPdN_data[i]);
  }
  free(p_data);
  free(prBC_data);
  free(vorticity_data);
  free(save_temp_data);
  free(new_order_data);
  free(rho_data);
  free(mu_data);
}

void INSData::init() {
  // DG_NP tmps
  // Advection
  F[0] = tmp_dg_np[0];
  F[1] = tmp_dg_np[1];
  F[2] = tmp_dg_np[2];
  F[3] = tmp_dg_np[3];
  // Pressure
  divVelT = tmp_dg_np[0];
  curlVel = tmp_dg_np[1];
  gradCurlVel[0] = tmp_dg_np[2];
  gradCurlVel[1] = tmp_dg_np[3];
  pRHS = tmp_dg_np[1];
  dpdx = tmp_dg_np[2];
  dpdy = tmp_dg_np[3];
  // Viscosity
  visRHS[0] = tmp_dg_np[0];
  visRHS[1] = tmp_dg_np[1];
  visTmp[0] = tmp_dg_np[2];
  visTmp[1] = tmp_dg_np[3];

  // DG_G_NP tmps
  // Advection
  gQ[0]   = tmp_dg_g_np[0];
  gQ[1]   = tmp_dg_g_np[1];
  flux[0] = tmp_dg_g_np[2];
  flux[1] = tmp_dg_g_np[3];
  // Pressure
  gP     = tmp_dg_g_np[0];
  pFluxX = tmp_dg_g_np[1];
  pFluxY = tmp_dg_g_np[2];
  gN[0]  = tmp_dg_g_np[0];
  gN[1]  = tmp_dg_g_np[1];
  gGradCurl[0] = tmp_dg_g_np[2];
  gGradCurl[1] = tmp_dg_g_np[3];
  gRho = tmp_dg_g_np[4];
  // Viscosity
  visBC[0] = tmp_dg_g_np[0];
  visBC[1] = tmp_dg_g_np[1];
}
