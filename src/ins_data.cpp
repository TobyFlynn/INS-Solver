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

  for(int i = 0; i < 10; i++) {
    tmp_dg_np_data[i] = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  }
  for(int i = 0; i < 4; i++) {
    tmp_dg_g_np_data[i] = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  }
  for(int i = 0; i < 2; i++) {
    Q_data[0][i]        = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    Q_data[1][i]        = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    QT_data[i]          = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    QTT_data[i]         = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    N_data[0][i]        = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    N_data[1][i]        = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    surf_ten_data[0][i] = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    surf_ten_data[1][i] = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
    dPdN_data[i]        = (double *)calloc(3 * DG_NPF * mesh->numCells, sizeof(double));
    tmp_dg_npf_data[i]  = (double *)calloc(3 * DG_NPF * mesh->numCells, sizeof(double));
  }
  p_data         = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  vorticity_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  save_temp_data = (double *)calloc(DG_SUB_CELLS * mesh->numCells, sizeof(double));
  nu_data        = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  gNu_data       = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  rho_data       = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  Dx_data    = (double *)calloc(DG_CUB_NP * DG_NP * mesh->numCells, sizeof(double));
  Dy_data    = (double *)calloc(DG_CUB_NP * DG_NP * mesh->numCells, sizeof(double));

  grx_data     = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  gsx_data     = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  gry_data     = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  gsy_data     = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  reverse_data = (int *)calloc(3 * mesh->numCells, sizeof(int));
  for(int i = 0; i < 3; i++) {
    mDx_data[i]  = (double *)calloc(DG_GF_NP * DG_NP * mesh->numCells, sizeof(double));
    mDy_data[i]  = (double *)calloc(DG_GF_NP * DG_NP * mesh->numCells, sizeof(double));
  }

  mDL_data  = (double *)calloc(DG_GF_NP * DG_NP * mesh->numEdges, sizeof(double));
  mDR_data  = (double *)calloc(DG_GF_NP * DG_NP * mesh->numEdges, sizeof(double));
  mDBC_data = (double *)calloc(DG_GF_NP * DG_NP * mesh->numBoundaryEdges, sizeof(double));
  pDL_data  = (double *)calloc(DG_GF_NP * DG_NP * mesh->numEdges, sizeof(double));
  pDR_data  = (double *)calloc(DG_GF_NP * DG_NP * mesh->numEdges, sizeof(double));
  gVPL_data = (double *)calloc(DG_GF_NP * DG_NP * mesh->numEdges, sizeof(double));
  gVPR_data = (double *)calloc(DG_GF_NP * DG_NP * mesh->numEdges, sizeof(double));

  for(int i = 0; i < 10; i++) {
    string name  = "tmp_dg_np" + to_string(i);
    tmp_dg_np[i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_dg_np_data[i], name.c_str());
  }
  for(int i = 0; i < 4; i++) {
    string name  = "tmp_dg_g_np" + to_string(i);
    tmp_dg_g_np[i] = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_dg_g_np_data[i], name.c_str());
  }
  for(int i = 0; i < 2; i++) {
    string name    = "Q0" + to_string(i);
    Q[0][i]        = op_decl_dat(mesh->cells, DG_NP, "double", Q_data[0][i], name.c_str());
    name           = "Q1" + to_string(i);
    Q[1][i]        = op_decl_dat(mesh->cells, DG_NP, "double", Q_data[1][i], name.c_str());
    name           = "QT" + to_string(i);
    QT[i]          = op_decl_dat(mesh->cells, DG_NP, "double", QT_data[i], name.c_str());
    name           = "QTT" + to_string(i);
    QTT[i]         = op_decl_dat(mesh->cells, DG_NP, "double", QTT_data[i], name.c_str());
    name           = "N0" + to_string(i);
    N[0][i]        = op_decl_dat(mesh->cells, DG_NP, "double", N_data[0][i], name.c_str());
    name           = "N1" + to_string(i);
    N[1][i]        = op_decl_dat(mesh->cells, DG_NP, "double", N_data[1][i], name.c_str());
    name           = "surf_ten0" + to_string(i);
    surf_ten[0][i] = op_decl_dat(mesh->cells, DG_NP, "double", surf_ten_data[0][i], name.c_str());
    name           = "surf_ten1" + to_string(i);
    surf_ten[1][i] = op_decl_dat(mesh->cells, DG_NP, "double", surf_ten_data[1][i], name.c_str());
    name           = "dPdN" + to_string(i);
    dPdN[i]        = op_decl_dat(mesh->cells, 3 * DG_NPF, "double", dPdN_data[i], name.c_str());
    name           = "tmp_dg_npf" + to_string(i);
    tmp_dg_npf[i]  = op_decl_dat(mesh->cells, 3 * DG_NPF, "double", tmp_dg_npf_data[i], name.c_str());
  }
  p         = op_decl_dat(mesh->cells, DG_NP, "double", p_data, "p");
  vorticity = op_decl_dat(mesh->cells, DG_NP, "double", vorticity_data, "vorticity");
  save_temp = op_decl_dat(mesh->cells, DG_SUB_CELLS, "double", save_temp_data, "save_temp");
  nu        = op_decl_dat(mesh->cells, DG_NP, "double", nu_data, "nu");
  gNu       = op_decl_dat(mesh->cells, DG_G_NP, "double", gNu_data, "gNu");
  rho       = op_decl_dat(mesh->cells, DG_NP, "double", rho_data, "rho");

  Dx    = op_decl_dat(mesh->cells, DG_CUB_NP * DG_NP, "double", Dx_data, "cub-Dx");
  Dy    = op_decl_dat(mesh->cells, DG_CUB_NP * DG_NP, "double", Dy_data, "cub-Dy");

  grx     = op_decl_dat(mesh->cells, DG_G_NP, "double", grx_data, "gauss-grx");
  gsx     = op_decl_dat(mesh->cells, DG_G_NP, "double", gsx_data, "gauss-gsx");
  gry     = op_decl_dat(mesh->cells, DG_G_NP, "double", gry_data, "gauss-gry");
  gsy     = op_decl_dat(mesh->cells, DG_G_NP, "double", gsy_data, "gauss-gsy");
  reverse = op_decl_dat(mesh->cells, 3, "int", reverse_data, "gauss-reverse");

  for(int i = 0; i < 3; i++) {
    string name = "gauss-mDx" + to_string(i);
    mDx[i]      = op_decl_dat(mesh->cells, DG_GF_NP * DG_NP, "double", mDx_data[i], name.c_str());
    name        = "gauss-mDy" + to_string(i);
    mDy[i]      = op_decl_dat(mesh->cells, DG_GF_NP * DG_NP, "double", mDy_data[i], name.c_str());
  }

  mDL  = op_decl_dat(mesh->edges, DG_GF_NP * DG_NP, "double", mDL_data, "mDL");
  mDR  = op_decl_dat(mesh->edges, DG_GF_NP * DG_NP, "double", mDR_data, "mDR");
  mDBC = op_decl_dat(mesh->bedges, DG_GF_NP * DG_NP, "double", mDBC_data, "mDBC");
  pDL  = op_decl_dat(mesh->edges, DG_GF_NP * DG_NP, "double", pDL_data, "pDL");
  pDR  = op_decl_dat(mesh->edges, DG_GF_NP * DG_NP, "double", pDR_data, "pDR");
  gVPL = op_decl_dat(mesh->edges, DG_GF_NP * DG_NP, "double", gVPL_data, "gVPL");
  gVPR = op_decl_dat(mesh->edges, DG_GF_NP * DG_NP, "double", gVPR_data, "gVPR");

  op_decl_const(1, "double", &reynolds);
  op_decl_const(1, "double", &froude);
  op_decl_const(1, "double", &weber);
  op_decl_const(1, "double", &nu0);
  op_decl_const(1, "double", &nu1);
  op_decl_const(1, "double", &rho0);
  op_decl_const(1, "double", &rho1);
  op_decl_const(3 * DG_NPF, "int", FMASK);
  op_decl_const(1, "double", &ic_u);
  op_decl_const(1, "double", &ic_v);
  op_decl_const(DG_CUB_NP, "double", cubW_g);
  op_decl_const(DG_CUB_NP * DG_NP, "double", cubV_g);
  op_decl_const(DG_CUB_NP * DG_NP, "double", cubVDr_g);
  op_decl_const(DG_CUB_NP * DG_NP, "double", cubVDs_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF0Dr_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF0Ds_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF1Dr_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF1Ds_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF2Dr_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF2Ds_g);
  op_decl_const(DG_GF_NP, "double", gaussW_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gFInterp0_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gFInterp1_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gFInterp2_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF0DrR_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF0DsR_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF1DrR_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF1DsR_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF2DrR_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gF2DsR_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gFInterp0R_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gFInterp1R_g);
  op_decl_const(DG_GF_NP * DG_NP, "double", gFInterp2R_g);
}

INSData::~INSData() {
  for(int i = 0; i < 10; i++) {
    free(tmp_dg_np_data[i]);
  }
  for(int i = 0; i < 4; i++) {
    free(tmp_dg_g_np_data[i]);
  }
  for(int i = 0; i < 2; i++) {
    free(Q_data[0][i]);
    free(Q_data[1][i]);
    free(QT_data[i]);
    free(QTT_data[i]);
    free(N_data[0][i]);
    free(N_data[1][i]);
    free(surf_ten_data[0][i]);
    free(surf_ten_data[1][i]);
    free(dPdN_data[i]);
    free(tmp_dg_npf_data[i]);
  }
  free(p_data);
  free(vorticity_data);
  free(save_temp_data);
  free(nu_data);
  free(gNu_data);
  free(rho_data);

  free(Dx_data);
  free(Dy_data);

  free(grx_data);
  free(gsx_data);
  free(gry_data);
  free(gsy_data);
  free(reverse_data);
  for(int i = 0; i < 3; i++) {
    free(mDx_data[i]);
    free(mDy_data[i]);
  }

  free(mDL_data);
  free(mDR_data);
  free(mDBC_data);
  free(pDL_data);
  free(pDR_data);
  free(gVPL_data);
  free(gVPR_data);
}

void INSData::init() {
  // Set up dats that share a storage dat
  F[0] = tmp_dg_np[0];
  F[1] = tmp_dg_np[1];
  F[2] = tmp_dg_np[2];
  F[3] = tmp_dg_np[3];

  divVelT        = tmp_dg_np[0];
  curlVel        = tmp_dg_np[1];
  gradCurlVel[0] = tmp_dg_np[2];
  gradCurlVel[1] = tmp_dg_np[3];

  pRHS       = tmp_dg_np[1];
  dpdx       = tmp_dg_np[2];
  dpdy       = tmp_dg_np[3];
  visRHS[0]  = tmp_dg_np[0];
  visRHS[1]  = tmp_dg_np[1];
  visTemp[0] = tmp_dg_np[2];
  visTemp[1] = tmp_dg_np[3];

  pFluxX  = tmp_dg_npf[0];
  pFluxY  = tmp_dg_npf[1];

  gQ[0]    = tmp_dg_g_np[0];
  gQ[1]    = tmp_dg_g_np[1];
  flux[0]  = tmp_dg_g_np[2];
  flux[1]  = tmp_dg_g_np[3];
  prBC     = tmp_dg_g_np[0];
  visBC[0] = tmp_dg_g_np[1];
  visBC[1] = tmp_dg_g_np[2];

  // Regular grid point init
  op_par_loop(init_nu_rho, "init_nu_rho", mesh->cells,
              op_arg_dat(nu,  -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rho, -1, OP_ID, DG_NP, "double", OP_WRITE));

  op2_gemv(false, DG_G_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_G_NP, nu, 0.0, gNu);

  // Cubature grid point init (needed for Poisson solver)
  // Initialise geometric factors for calcuating grad matrix
  op2_gemv(false, DG_CUB_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::CUB_VDR), DG_CUB_NP, mesh->x, 0.0, mesh->cubature->op_tmp[0]);
  op2_gemv(false, DG_CUB_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::CUB_VDS), DG_CUB_NP, mesh->x, 0.0, mesh->cubature->op_tmp[1]);
  op2_gemv(false, DG_CUB_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::CUB_VDR), DG_CUB_NP, mesh->y, 0.0, mesh->cubature->op_tmp[2]);
  op2_gemv(false, DG_CUB_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::CUB_VDS), DG_CUB_NP, mesh->y, 0.0, mesh->cubature->op_tmp[3]);

  // The Dx and Dy dats contain matrices that are used when calculating the 1st term of Eqn. 10 in Karakus et al.
  op_par_loop(init_cubature_grad, "init_cubature_grad", mesh->cells,
              op_arg_dat(mesh->cubature->op_tmp[0], -1, OP_ID, DG_CUB_NP, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[1], -1, OP_ID, DG_CUB_NP, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[2], -1, OP_ID, DG_CUB_NP, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[3], -1, OP_ID, DG_CUB_NP, "double", OP_RW),
              op_arg_dat(Dx, -1, OP_ID, DG_CUB_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(Dy, -1, OP_ID, DG_CUB_NP * DG_NP, "double", OP_WRITE));
  // Dx and Dy are col-major at this point

  /*****************************************************************************
  *
  * Below contains code used to calculate matrices which are later used to
  * calculate terms 2, 3 and 4 of Eqn. 10 in Karakus et al. You can ignore most
  * of these dats, most are just temp dats, the dats holding the final matrices
  * are mDL, mDR, mDBC, pDL, pDR, gVPL and gVPR.
  *
  * Looking at the following code from the Hesthaven and Warburton textbook can
  * help to understand these matrices and how they are calculated:
  * https://github.com/tcew/nodal-dg/blob/master/Codes1.1/Codes2D/CurvedPoissonIPDG2D.m
  * mDL, mDR and mDBC correspond to gDnM in the MATLAB code. pDL and pDR to gDnP.
  * gVPL and gVPR to gVP. L dats belong to the cell to the left of the edge, R
  * to the right cell.
  *
  *****************************************************************************/

  // Gauss grid point init
  // Check which edges will require matrices to be 'reverse'
  op_par_loop(gauss_reverse, "gauss_reverse", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->nodeX,   -2, mesh->edge2cells, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY,   -2, mesh->edge2cells, 3, "double", OP_READ),
              op_arg_dat(reverse,       -2, mesh->edge2cells, 3, "int", OP_INC));

  // Calculate geometric factors used when constructing gradient matrices
  init_gauss_grad_blas(mesh, this);

  // [gDxM, gDyM] = PhysDmatrices2D(x(:,k1), y(:,k1),gVM);
  op_par_loop(init_gauss_grad, "init_gauss_grad", mesh->cells,
              op_arg_dat(grx,    -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(gsx,    -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(gry,    -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(gsy,    -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(mDx[0], -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDy[0], -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDx[1], -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDy[1], -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDx[2], -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDy[2], -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE));

  // gDnM = gnx*gDxM + gny*gDyM;
  op_par_loop(init_gauss_grad3, "init_gauss_grad3", mesh->edges,
              op_arg_dat(mesh->edgeNum,  -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 1, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 1, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mDx[0], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[0], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[0], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[0], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[1], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[1], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[1], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[1], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[2], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[2], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[2], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[2], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDL, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDR, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE));

  op_par_loop(init_gauss_grad4, "init_gauss_grad4", mesh->bedges,
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mDx[0], 0, mesh->bedge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[0], 0, mesh->bedge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[1], 0, mesh->bedge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[1], 0, mesh->bedge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[2], 0, mesh->bedge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[2], 0, mesh->bedge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDBC, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE));

  // Calculate geometric factors for grad matrices used by neighbours
  // Matrices are calculated locally, then copied to neighbour elements
  init_gauss_grad_neighbour_blas(mesh, this);

  // [gDxP, gDyP] = PhysDmatrices2D(x(:,k2), y(:,k2),gVP);
  op_par_loop(init_gauss_grad_neighbour, "init_gauss_grad_neighbour", mesh->cells,
              op_arg_dat(reverse, -1, OP_ID, 3, "int", OP_READ),
              op_arg_dat(grx,     -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(gsx,     -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(gry,     -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(gsy,     -1, OP_ID, DG_G_NP, "double", OP_RW),
              op_arg_dat(mDx[0],  -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDy[0],  -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDx[1],  -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDy[1],  -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDx[2],  -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(mDy[2],  -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE));

  // gDnP = gnx*gDxP + gny*gDyP;
  op_par_loop(init_gauss_grad5, "init_gauss_grad5", mesh->edges,
              op_arg_dat(mesh->edgeNum,  -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 1, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 1, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mDx[0], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[0], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[0], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[0], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[1], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[1], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[1], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[1], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[2], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDx[2], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[2], 0, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mDy[2], 1, mesh->edge2cells, DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(pDL, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(pDR, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE));

  op_par_loop(gauss_gfi_faces2, "gauss_gfi_faces2", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(gVPL, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(gVPR, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE));
}
