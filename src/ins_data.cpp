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

  double *tmp_np = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  double *tmp_g_np = (double *)calloc(DG_G_NP * mesh->cells->size, sizeof(double));
  double *tmp_sub_cell = (double *)calloc(DG_SUB_CELLS * mesh->cells->size, sizeof(double));
  int* tmp_int_1 = (int *)calloc(mesh->cells->size, sizeof(int));

  // Declare OP2 datasets
  for(int i = 0; i < 10; i++) {
    string name    = "tmp_dg_np" + to_string(i);
    tmp_dg_np[i]   = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
  }
  for(int i = 0; i < 5; i++) {
    string name    = "tmp_dg_g_np" + to_string(i);
    tmp_dg_g_np[i] = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_g_np, name.c_str());
  }
  for(int i = 0; i < 2; i++) {
    string name = "Q0" + to_string(i);
    Q[0][i]     = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name        = "Q1" + to_string(i);
    Q[1][i]     = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name        = "QT" + to_string(i);
    QT[i]       = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name        = "QTT" + to_string(i);
    QTT[i]      = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name        = "N0" + to_string(i);
    N[0][i]     = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name        = "N1" + to_string(i);
    N[1][i]     = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name        = "dPdN" + to_string(i);
    dPdN[i]     = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_g_np, name.c_str());
  }
  p         = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "p");
  prBC      = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_g_np, "prBC");
  vorticity = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "vorticity");
  save_temp = op_decl_dat(mesh->cells, DG_SUB_CELLS, "double", tmp_sub_cell, "save_temp");
  new_order = op_decl_dat(mesh->cells, 1, "int", tmp_int_1, "new_order");
  rho       = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "rho");
  mu        = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "mu");

  free(tmp_int_1);
  free(tmp_sub_cell);
  free(tmp_g_np);
  free(tmp_np);

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
