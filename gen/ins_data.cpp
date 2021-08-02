#include "ins_data.h"

#include "op_seq.h"

#include <string>
#include <memory>

#include "dg_blas_calls.h"
#include "dg_compiler_defs.h"
#include "blas_calls.h"

#include "constants.h"

using namespace std;

INSData::INSData(DGMesh *m) {
  mesh = m;

  for(int i = 0; i < 4; i++) {
    tmp_dg_np_data[i] = (double *)calloc(10 * mesh->numCells, sizeof(double));
  }
  for(int i = 0; i < 2; i++) {
    Q_data[0][i]       = (double *)calloc(10 * mesh->numCells, sizeof(double));
    Q_data[1][i]       = (double *)calloc(10 * mesh->numCells, sizeof(double));
    QT_data[i]         = (double *)calloc(10 * mesh->numCells, sizeof(double));
    QTT_data[i]        = (double *)calloc(10 * mesh->numCells, sizeof(double));
    N_data[0][i]       = (double *)calloc(10 * mesh->numCells, sizeof(double));
    N_data[1][i]       = (double *)calloc(10 * mesh->numCells, sizeof(double));
    dPdN_data[i]       = (double *)calloc(3 * 4 * mesh->numCells, sizeof(double));
    visBC_data[i]      = (double *)calloc(18 * mesh->numCells, sizeof(double));
    tmp_dg_npf_data[i] = (double *)calloc(3 * 4 * mesh->numCells, sizeof(double));
  }
  p_data         = (double *)calloc(10 * mesh->numCells, sizeof(double));
  prBC_data      = (double *)calloc(18 * mesh->numCells, sizeof(double));
  vorticity_data = (double *)calloc(10 * mesh->numCells, sizeof(double));
  save_temp_data = (double *)calloc(9 * mesh->numCells, sizeof(double));
  nu_data        = (double *)calloc(10 * mesh->numCells, sizeof(double));
  gNu_data       = (double *)calloc(18 * mesh->numCells, sizeof(double));
  rho_data       = (double *)calloc(10 * mesh->numCells, sizeof(double));
  pFluxX_data    = (double *)calloc(3 * 4 * mesh->numCells, sizeof(double));
  pFluxY_data    = (double *)calloc(3 * 4 * mesh->numCells, sizeof(double));

  Dx_data    = (double *)calloc(36 * 10 * mesh->numCells, sizeof(double));
  Dy_data    = (double *)calloc(36 * 10 * mesh->numCells, sizeof(double));
  cOP_data   = (double *)calloc(10 * 10 * mesh->numCells, sizeof(double));
  temp_data  = (double *)calloc(36 * 10 * mesh->numCells, sizeof(double));
  temp2_data = (double *)calloc(36 * 10 * mesh->numCells, sizeof(double));

  grx_data     = (double *)calloc(18 * mesh->numCells, sizeof(double));
  gsx_data     = (double *)calloc(18 * mesh->numCells, sizeof(double));
  gry_data     = (double *)calloc(18 * mesh->numCells, sizeof(double));
  gsy_data     = (double *)calloc(18 * mesh->numCells, sizeof(double));
  tau_data     = (double *)calloc(3 * mesh->numCells, sizeof(double));
  reverse_data = (int *)calloc(3 * mesh->numCells, sizeof(int));
  for(int i = 0; i < 3; i++) {
    mDx_data[i]  = (double *)calloc(6 * 10 * mesh->numCells, sizeof(double));
    mDy_data[i]  = (double *)calloc(6 * 10 * mesh->numCells, sizeof(double));
    pDx_data[i]  = (double *)calloc(6 * 10 * mesh->numCells, sizeof(double));
    pDy_data[i]  = (double *)calloc(6 * 10 * mesh->numCells, sizeof(double));
    mD_data[i]   = (double *)calloc(6 * 10 * mesh->numCells, sizeof(double));
    pD_data[i]   = (double *)calloc(6 * 10 * mesh->numCells, sizeof(double));
    gOP_data[i]  = (double *)calloc(10 * 10 * mesh->numCells, sizeof(double));
    gOPf_data[i] = (double *)calloc(10 * 10 * mesh->numCells, sizeof(double));
  }

  for(int i = 0; i < 4; i++) {
    string name  = "tmp_dg_np" + to_string(i);
    tmp_dg_np[i] = op_decl_dat(mesh->cells, 10, "double", tmp_dg_np_data[i], name.c_str());
  }
  for(int i = 0; i < 2; i++) {
    string name   = "Q0" + to_string(i);
    Q[0][i]       = op_decl_dat(mesh->cells, 10, "double", Q_data[0][i], name.c_str());
    name          = "Q1" + to_string(i);
    Q[1][i]       = op_decl_dat(mesh->cells, 10, "double", Q_data[1][i], name.c_str());
    name          = "QT" + to_string(i);
    QT[i]         = op_decl_dat(mesh->cells, 10, "double", QT_data[i], name.c_str());
    name          = "QTT" + to_string(i);
    QTT[i]        = op_decl_dat(mesh->cells, 10, "double", QTT_data[i], name.c_str());
    name          = "N0" + to_string(i);
    N[0][i]       = op_decl_dat(mesh->cells, 10, "double", N_data[0][i], name.c_str());
    name          = "N1" + to_string(i);
    N[1][i]       = op_decl_dat(mesh->cells, 10, "double", N_data[1][i], name.c_str());
    name          = "dPdN" + to_string(i);
    dPdN[i]       = op_decl_dat(mesh->cells, 3 * 4, "double", dPdN_data[i], name.c_str());
    name          = "visBC" + to_string(i);
    visBC[i]      = op_decl_dat(mesh->cells, 18, "double", visBC_data[i], name.c_str());
    name          = "tmp_dg_npf" + to_string(i);
    tmp_dg_npf[i] = op_decl_dat(mesh->cells, 3 * 4, "double", tmp_dg_npf_data[i], name.c_str());
  }
  p         = op_decl_dat(mesh->cells, 10, "double", p_data, "p");
  prBC      = op_decl_dat(mesh->cells, 18, "double", prBC_data, "prBC");
  vorticity = op_decl_dat(mesh->cells, 10, "double", vorticity_data, "vorticity");
  save_temp = op_decl_dat(mesh->cells, 9, "double", save_temp_data, "save_temp");
  nu        = op_decl_dat(mesh->cells, 10, "double", nu_data, "nu");
  gNu       = op_decl_dat(mesh->cells, 18, "double", gNu_data, "gNu");
  rho       = op_decl_dat(mesh->cells, 10, "double", rho_data, "rho");
  pFluxX    = op_decl_dat(mesh->cells, 3 * 4, "double", pFluxX_data, "pX");
  pFluxY    = op_decl_dat(mesh->cells, 3 * 4, "double", pFluxY_data, "pY");

  Dx    = op_decl_dat(mesh->cells, 36 * 10, "double", Dx_data, "cub-Dx");
  Dy    = op_decl_dat(mesh->cells, 36 * 10, "double", Dy_data, "cub-Dy");
  cOP   = op_decl_dat(mesh->cells, 10 * 10, "double", cOP_data, "cub-OP");
  temp  = op_decl_dat(mesh->cells, 36 * 10, "double", temp_data, "cub-temp");
  temp2 = op_decl_dat(mesh->cells, 36 * 10, "double", temp2_data, "cub-temp2");

  grx     = op_decl_dat(mesh->cells, 18, "double", grx_data, "gauss-grx");
  gsx     = op_decl_dat(mesh->cells, 18, "double", gsx_data, "gauss-gsx");
  gry     = op_decl_dat(mesh->cells, 18, "double", gry_data, "gauss-gry");
  gsy     = op_decl_dat(mesh->cells, 18, "double", gsy_data, "gauss-gsy");
  tau     = op_decl_dat(mesh->cells, 3, "double", tau_data, "gauss-tau");
  reverse = op_decl_dat(mesh->cells, 3, "int", reverse_data, "gauss-reverse");

  for(int i = 0; i < 3; i++) {
    string name = "gauss-mDx" + to_string(i);
    mDx[i]      = op_decl_dat(mesh->cells, 6 * 10, "double", mDx_data[i], name.c_str());
    name        = "gauss-mDy" + to_string(i);
    mDy[i]      = op_decl_dat(mesh->cells, 6 * 10, "double", mDy_data[i], name.c_str());
    name        = "gauss-pDx" + to_string(i);
    pDx[i]      = op_decl_dat(mesh->cells, 6 * 10, "double", pDx_data[i], name.c_str());
    name        = "gauss-pDy" + to_string(i);
    pDy[i]      = op_decl_dat(mesh->cells, 6 * 10, "double", pDy_data[i], name.c_str());
    name        = "gauss-mD" + to_string(i);
    mD[i]       = op_decl_dat(mesh->cells, 6 * 10, "double", mD_data[i], name.c_str());
    name        = "gauss-pD" + to_string(i);
    pD[i]       = op_decl_dat(mesh->cells, 6 * 10, "double", pD_data[i], name.c_str());
    name        = "gauss-OP" + to_string(i);
    gOP[i]      = op_decl_dat(mesh->cells, 10 * 10, "double", gOP_data[i], name.c_str());
    name        = "gauss-OPf" + to_string(i);
    gOPf[i]     = op_decl_dat(mesh->cells, 10 * 10, "double", gOPf_data[i], name.c_str());
  }

  op_decl_const(1, "double", &reynolds);
  op_decl_const(1, "double", &froude);
  op_decl_const(1, "double", &weber);
  op_decl_const(1, "double", &nu0);
  op_decl_const(1, "double", &nu1);
  op_decl_const(1, "double", &rho0);
  op_decl_const(1, "double", &rho1);
  op_decl_const(3 * 4, "int", FMASK);
  op_decl_const(1, "double", &ic_u);
  op_decl_const(1, "double", &ic_v);
  op_decl_const(36, "double", cubW_g);
  op_decl_const(36 * 10, "double", cubV_g);
  op_decl_const(36 * 10, "double", cubVDr_g);
  op_decl_const(36 * 10, "double", cubVDs_g);
  op_decl_const(6 * 10, "double", gF0Dr_g);
  op_decl_const(6 * 10, "double", gF0Ds_g);
  op_decl_const(6 * 10, "double", gF1Dr_g);
  op_decl_const(6 * 10, "double", gF1Ds_g);
  op_decl_const(6 * 10, "double", gF2Dr_g);
  op_decl_const(6 * 10, "double", gF2Ds_g);
  op_decl_const(6, "double", gaussW_g);
  op_decl_const(6 * 10, "double", gFInterp0_g);
  op_decl_const(6 * 10, "double", gFInterp1_g);
  op_decl_const(6 * 10, "double", gFInterp2_g);
  op_decl_const(6 * 10, "double", gF0DrR_g);
  op_decl_const(6 * 10, "double", gF0DsR_g);
  op_decl_const(6 * 10, "double", gF1DrR_g);
  op_decl_const(6 * 10, "double", gF1DsR_g);
  op_decl_const(6 * 10, "double", gF2DrR_g);
  op_decl_const(6 * 10, "double", gF2DsR_g);
  op_decl_const(6 * 10, "double", gFInterp0R_g);
  op_decl_const(6 * 10, "double", gFInterp1R_g);
  op_decl_const(6 * 10, "double", gFInterp2R_g);
}

INSData::~INSData() {
  for(int i = 0; i < 4; i++) {
    free(tmp_dg_np_data[i]);
  }
  for(int i = 0; i < 2; i++) {
    free(Q_data[0][i]);
    free(Q_data[1][i]);
    free(QT_data[i]);
    free(QTT_data[i]);

    free(N_data[0][i]);
    free(N_data[1][i]);
    free(dPdN_data[i]);
    free(visBC_data[i]);
    free(tmp_dg_npf_data[i]);
  }
  free(p_data);
  free(prBC_data);
  free(vorticity_data);
  free(save_temp_data);
  free(nu_data);
  free(gNu_data);
  free(rho_data);
  free(pFluxX_data);
  free(pFluxY_data);

  free(Dx_data);
  free(Dy_data);
  free(cOP_data);
  free(temp_data);
  free(temp2_data);

  free(grx_data);
  free(gsx_data);
  free(gry_data);
  free(gsy_data);
  free(tau_data);
  free(reverse_data);
  for(int i = 0; i < 3; i++) {
    free(mDx_data[i]);
    free(mDy_data[i]);
    free(pDx_data[i]);
    free(pDy_data[i]);
    free(mD_data[i]);
    free(pD_data[i]);
    free(gOP_data[i]);
    free(gOPf_data[i]);
  }
}

void INSData::init() {
  // Set up dats that share a storage dat
  F[0] = tmp_dg_np[0]; F[1] = tmp_dg_np[1];
  F[2] = tmp_dg_np[2]; F[3] = tmp_dg_np[3];
  divVelT = tmp_dg_np[0];
  curlVel = tmp_dg_np[1];
  gradCurlVel[0] = tmp_dg_np[2];
  gradCurlVel[1] = tmp_dg_np[3];
  pRHS = tmp_dg_np[1];
  dpdx = tmp_dg_np[2];
  dpdy = tmp_dg_np[3];
  visRHS[0] = tmp_dg_np[0];
  visRHS[1] = tmp_dg_np[1];
  flux[0] = tmp_dg_npf[0];
  flux[1] = tmp_dg_npf[1];
  pFluxX = tmp_dg_npf[0];
  pFluxY = tmp_dg_npf[1];

  // Regular grid point init
  op_par_loop(init_nu_rho, "init_nu_rho", mesh->cells,
              op_arg_dat(nu,  -1, OP_ID, 10, "double", OP_WRITE),
              op_arg_dat(rho, -1, OP_ID, 10, "double", OP_WRITE));

  op2_gemv(true, 18, 10, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), 10, nu, 0.0, gNu);

  // Cubature grid point init (needed for Poisson solver)
  // Initialise geometric factors for calcuating grad matrix
  op2_gemv(true, 36, 10, 1.0, constants->get_ptr(DGConstants::CUB_VDR), 10, mesh->x, 0.0, mesh->cubature->op_tmp[0]);
  op2_gemv(true, 36, 10, 1.0, constants->get_ptr(DGConstants::CUB_VDS), 10, mesh->x, 0.0, mesh->cubature->op_tmp[1]);
  op2_gemv(true, 36, 10, 1.0, constants->get_ptr(DGConstants::CUB_VDR), 10, mesh->y, 0.0, mesh->cubature->op_tmp[2]);
  op2_gemv(true, 36, 10, 1.0, constants->get_ptr(DGConstants::CUB_VDS), 10, mesh->y, 0.0, mesh->cubature->op_tmp[3]);

  op_par_loop(init_cubature_grad, "init_cubature_grad", mesh->cells,
              op_arg_dat(mesh->cubature->op_tmp[0], -1, OP_ID, 36, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[1], -1, OP_ID, 36, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[2], -1, OP_ID, 36, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[3], -1, OP_ID, 36, "double", OP_RW),
              op_arg_dat(Dx, -1, OP_ID, 36 * 10, "double", OP_WRITE),
              op_arg_dat(Dy, -1, OP_ID, 36 * 10, "double", OP_WRITE));
  // Dx and Dy are row-major at this point

  // Calculate Cubature OP (contribution of Cubature points to Poisson matrix)
  op_par_loop(init_cubature_OP, "init_cubature_OP", mesh->cells,
              op_arg_dat(mesh->cubature->J, -1, OP_ID, 36, "double", OP_READ),
              op_arg_dat(Dx,    -1, OP_ID, 36 * 10, "double", OP_READ),
              op_arg_dat(Dy,    -1, OP_ID, 36 * 10, "double", OP_READ),
              op_arg_dat(temp,  -1, OP_ID, 36 * 10, "double", OP_WRITE),
              op_arg_dat(temp2, -1, OP_ID, 36 * 10, "double", OP_WRITE));
  // Temp and temp2 are in row-major at this point
  op2_gemm_batch(false, true, 10, 10, 36, 1.0, Dx, 10, temp, 10, 0.0, cOP, 10);
  op2_gemm_batch(false, true, 10, 10, 36, 1.0, Dy, 10, temp2, 10, 1.0, cOP, 10);
  // OP is in col-major at this point

  // Gauss grid point init
  // Check which edges will require matrices to be 'reverse'
  op_par_loop(gauss_reverse, "gauss_reverse", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->nodeX,   -2, mesh->edge2cells, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY,   -2, mesh->edge2cells, 3, "double", OP_READ),
              op_arg_dat(reverse,       -2, mesh->edge2cells, 3, "int", OP_INC));

  // Calculate tau (used when constructing the Poisson matrix)
  op_par_loop(gauss_tau, "gauss_tau", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fscale,  -2, mesh->edge2cells, 3 * 4, "double", OP_READ),
              op_arg_dat(tau,           -2, mesh->edge2cells, 3, "double", OP_INC));

  op_par_loop(gauss_tau_bc, "gauss_tau_bc", mesh->bedges,
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->fscale,    0, mesh->bedge2cells, 3 * 4, "double", OP_READ),
              op_arg_dat(tau,             0, mesh->bedge2cells, 3, "double", OP_INC));

  // Calculate geometric factors used when constructing gradient matrices
  init_gauss_grad_blas(mesh, this);

  op_par_loop(init_gauss_grad, "init_gauss_grad", mesh->cells,
              op_arg_dat(grx,    -1, OP_ID, 18, "double", OP_RW),
              op_arg_dat(gsx,    -1, OP_ID, 18, "double", OP_RW),
              op_arg_dat(gry,    -1, OP_ID, 18, "double", OP_RW),
              op_arg_dat(gsy,    -1, OP_ID, 18, "double", OP_RW),
              op_arg_dat(mDx[0], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDy[0], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDx[1], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDy[1], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDx[2], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDy[2], -1, OP_ID, 6 * 10, "double", OP_WRITE));

  // Construct gradient matrices (1 for each face)
  op_par_loop(init_gauss_grad2, "init_gauss_grad2", mesh->cells,
              op_arg_dat(mesh->gauss->nx, -1, OP_ID, 18, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -1, OP_ID, 18, "double", OP_READ),
              op_arg_dat(mDx[0], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(mDy[0], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(mDx[1], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(mDy[1], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(mDx[2], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(mDy[2], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(mD[0],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mD[1],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mD[2],  -1, OP_ID, 6 * 10, "double", OP_WRITE));

  // Calculate geometric factors for grad matrices used by neighbours
  // Matrices are calculated locally, then copied to neighbour elements
  init_gauss_grad_neighbour_blas(mesh, this);

  op_par_loop(init_gauss_grad_neighbour, "init_gauss_grad_neighbour", mesh->cells,
              op_arg_dat(reverse, -1, OP_ID, 3, "int", OP_READ),
              op_arg_dat(grx,     -1, OP_ID, 18, "double", OP_RW),
              op_arg_dat(gsx,     -1, OP_ID, 18, "double", OP_RW),
              op_arg_dat(gry,     -1, OP_ID, 18, "double", OP_RW),
              op_arg_dat(gsy,     -1, OP_ID, 18, "double", OP_RW),
              op_arg_dat(mDx[0],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDy[0],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDx[1],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDy[1],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDx[2],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDy[2],  -1, OP_ID, 6 * 10, "double", OP_WRITE));

  // Copy x and y grad matrices to neighbours
  op_par_loop(gauss_grad_faces, "gauss_grad_faces", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mDx[0], -2, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(mDy[0], -2, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(mDx[1], -2, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(mDy[1], -2, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(mDx[2], -2, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(mDy[2], -2, mesh->edge2cells, 6 * 10, "double", OP_READ),
              op_arg_dat(pDx[0], -2, mesh->edge2cells, 6 * 10, "double", OP_INC),
              op_arg_dat(pDy[0], -2, mesh->edge2cells, 6 * 10, "double", OP_INC),
              op_arg_dat(pDx[1], -2, mesh->edge2cells, 6 * 10, "double", OP_INC),
              op_arg_dat(pDy[1], -2, mesh->edge2cells, 6 * 10, "double", OP_INC),
              op_arg_dat(pDx[2], -2, mesh->edge2cells, 6 * 10, "double", OP_INC),
              op_arg_dat(pDy[2], -2, mesh->edge2cells, 6 * 10, "double", OP_INC));

  // Calculate final neighbour grad matrix
  op_par_loop(init_gauss_grad2, "init_gauss_grad2", mesh->cells,
              op_arg_dat(mesh->gauss->nx, -1, OP_ID, 18, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -1, OP_ID, 18, "double", OP_READ),
              op_arg_dat(pDx[0], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(pDy[0], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(pDx[1], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(pDy[1], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(pDx[2], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(pDy[2], -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(pD[0],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(pD[1],  -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(pD[2],  -1, OP_ID, 6 * 10, "double", OP_WRITE));

  // Calculate Gauss OP for each face (local contribution of face in Poisson matrix)
  // Face 0 temps: mDx, Face 1 temps: mDy, Face 2 temps: pDx
  op_par_loop(gauss_op, "gauss_op", mesh->cells,
              op_arg_dat(tau, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -1, OP_ID, 18, "double", OP_READ),
              // Face 0
              op_arg_dat(mD[0],  -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(mDx[0], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDx[1], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDx[2], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              // Face 1
              op_arg_dat(mD[1],  -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(mDy[0], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDy[1], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(mDy[2], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              // Face 2
              op_arg_dat(mD[2],  -1, OP_ID, 6 * 10, "double", OP_READ),
              op_arg_dat(pDx[0], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(pDx[1], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(pDx[2], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              // Reset dats for OPf
              op_arg_dat(pDy[0], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(pDy[1], -1, OP_ID, 6 * 10, "double", OP_WRITE),
              op_arg_dat(pDy[2], -1, OP_ID, 6 * 10, "double", OP_WRITE));

  op2_gemm(true, true, 10, 10, 6, 1.0, mDx[0], 6, constants->get_ptr(DGConstants::GAUSS_FINTERP0), 10, 0.0, gOP[0], 10);
  op2_gemm_batch(true, true, 10, 10, 6, -1.0, mDx[1], 6, mD[0], 15, 1.0, gOP[0], 10);
  op2_gemm(true, true, 10, 10, 6, -1.0, mDx[2], 6, constants->get_ptr(DGConstants::GAUSS_FINTERP0), 10, 1.0, gOP[0], 10);

  op2_gemm(true, true, 10, 10, 6, 1.0, mDy[0], 6, constants->get_ptr(DGConstants::GAUSS_FINTERP1), 10, 0.0, gOP[1], 10);
  op2_gemm_batch(true, true, 10, 10, 6, -1.0, mDy[1], 6, mD[1], 10, 1.0, gOP[1], 10);
  op2_gemm(true, true, 10, 10, 6, -1.0, mDy[2], 6, constants->get_ptr(DGConstants::GAUSS_FINTERP1), 10, 1.0, gOP[1], 10);

  op2_gemm(true, true, 10, 10, 6, 1.0, pDx[0], 6, constants->get_ptr(DGConstants::GAUSS_FINTERP2), 10, 0.0, gOP[2], 10);
  op2_gemm_batch(true, true, 10, 10, 6, -1.0, pDx[1], 6, mD[2], 10, 1.0, gOP[2], 10);
  op2_gemm(true, true, 10, 10, 6, -1.0, pDx[2], 6, constants->get_ptr(DGConstants::GAUSS_FINTERP2), 10, 1.0, gOP[2], 10);

  // Calculate Gauss OPf for each face (contribution to neighbouring element in Poisson matrix)
  op_par_loop(gauss_gfi_faces, "gauss_gfi_faces", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(pDy[0], -2, mesh->edge2cells, 6 * 10, "double", OP_INC),
              op_arg_dat(pDy[1], -2, mesh->edge2cells, 6 * 10, "double", OP_INC),
              op_arg_dat(pDy[2], -2, mesh->edge2cells, 6 * 10, "double", OP_INC));

  op2_gemm_batch(true, true, 10, 10, 6, 1.0, mDx[0], 6, pDy[0], 10, 0.0, gOPf[0], 10);
  op2_gemm_batch(true, true, 10, 10, 6, 1.0, mDx[1], 6, pD[0], 10, 1.0, gOPf[0], 10);
  op2_gemm_batch(true, true, 10, 10, 6, -1.0, mDx[2], 6, pDy[0], 10, 1.0, gOPf[0], 10);

  op2_gemm_batch(true, true, 10, 10, 6, 1.0, mDy[0], 6, pDy[1], 10, 0.0, gOPf[1], 10);
  op2_gemm_batch(true, true, 10, 10, 6, 1.0, mDy[1], 6, pD[1], 10, 1.0, gOPf[1], 10);
  op2_gemm_batch(true, true, 10, 10, 6, -1.0, mDy[2], 6, pDy[1], 10, 1.0, gOPf[1], 10);

  op2_gemm_batch(true, true, 10, 10, 6, 1.0, pDx[0], 6, pDy[2], 10, 0.0, gOPf[2], 10);
  op2_gemm_batch(true, true, 10, 10, 6, 1.0, pDx[1], 6, pD[2], 10, 1.0, gOPf[2], 10);
  op2_gemm_batch(true, true, 10, 10, 6, -1.0, pDx[2], 6, pDy[2], 10, 1.0, gOPf[2], 10);

  // Applying the correct factors to OP and OPf is done when constructing the Poisson matrix
}
