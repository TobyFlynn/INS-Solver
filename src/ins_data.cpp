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
  // Initialise memory
  for(int i = 0; i < 4; i++) {
    F_data[i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    div_data[i] = (double *)calloc(15 * mesh->numCells, sizeof(double));
  }
  for(int i = 0; i < 2; i++) {
    Q_data[0][i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    Q_data[1][i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    QT_data[i]     = (double *)calloc(15 * mesh->numCells, sizeof(double));
    QTT_data[i]    = (double *)calloc(15 * mesh->numCells, sizeof(double));
    N_data[0][i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    N_data[1][i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    exQ_data[i]    = (double *)calloc(15 * mesh->numCells, sizeof(double));
    flux_data[i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    dPdN_data[i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    visRHS_data[i] = (double *)calloc(15 * mesh->numCells, sizeof(double));
    visBC_data[i]  = (double *)calloc(21 * mesh->numCells, sizeof(double));
    dQdx_data[i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    dQdy_data[i]   = (double *)calloc(15 * mesh->numCells, sizeof(double));
    gradCurlVel_data[i] = (double *)calloc(15 * mesh->numCells, sizeof(double));
  }
  divVelT_data   = (double *)calloc(15 * mesh->numCells, sizeof(double));
  curlVel_data   = (double *)calloc(15 * mesh->numCells, sizeof(double));
  pRHS_data      = (double *)calloc(15 * mesh->numCells, sizeof(double));
  pRHSex_data    = (double *)calloc(15 * mesh->numCells, sizeof(double));
  p_data         = (double *)calloc(15 * mesh->numCells, sizeof(double));
  dpdx_data      = (double *)calloc(15 * mesh->numCells, sizeof(double));
  dpdy_data      = (double *)calloc(15 * mesh->numCells, sizeof(double));
  prBC_data      = (double *)calloc(21 * mesh->numCells, sizeof(double));
  vorticity_data = (double *)calloc(15 * mesh->numCells, sizeof(double));
  save_temp_data = (double *)calloc(16 * mesh->numCells, sizeof(double));

  // Declare OP2 datasets
  for(int i = 0; i < 4; i++) {
    string Fname = "F" + to_string(i);
    F[i] = op_decl_dat(mesh->cells, 15, "double", F_data[i], Fname.c_str());
    string divname = "div" + to_string(i);
    div[i] = op_decl_dat(mesh->cells, 15, "double", div_data[i], divname.c_str());
  }
  for(int i = 0; i < 2; i++) {
    string Qname = "Q0" + to_string(i);
    Q[0][i] = op_decl_dat(mesh->cells, 15, "double", Q_data[0][i], Qname.c_str());
    Qname = "Q1" + to_string(i);
    Q[1][i] = op_decl_dat(mesh->cells, 15, "double", Q_data[1][i], Qname.c_str());
    Qname = "QT" + to_string(i);
    QT[i] = op_decl_dat(mesh->cells, 15, "double", QT_data[i], Qname.c_str());
    Qname = "QTT" + to_string(i);
    QTT[i] = op_decl_dat(mesh->cells, 15, "double", QTT_data[i], Qname.c_str());
    string Nname = "N0" + to_string(i);
    N[0][i] = op_decl_dat(mesh->cells, 15, "double", N_data[0][i], Nname.c_str());
    Nname = "N1" + to_string(i);
    N[1][i] = op_decl_dat(mesh->cells, 15, "double", N_data[1][i], Nname.c_str());
    string exQname = "exQ" + to_string(i);
    exQ[i] = op_decl_dat(mesh->cells, 15, "double", exQ_data[i], exQname.c_str());
    string fluxname = "flux" + to_string(i);
    flux[i] = op_decl_dat(mesh->cells, 15, "double", flux_data[i], fluxname.c_str());
    string gradCurlVelname = "gradCurlVel" + to_string(i);
    gradCurlVel[i] = op_decl_dat(mesh->cells, 15, "double", gradCurlVel_data[i], gradCurlVelname.c_str());
    string dPdNname = "dPdN" + to_string(i);
    dPdN[i] = op_decl_dat(mesh->cells, 15, "double", dPdN_data[i], dPdNname.c_str());
    string visRHSname = "visRHS" + to_string(i);
    visRHS[i] = op_decl_dat(mesh->cells, 15, "double", visRHS_data[i], visRHSname.c_str());
    string dQdxname = "dQdx" + to_string(i);
    dQdx[i] = op_decl_dat(mesh->cells, 15, "double", dQdx_data[i], dQdxname.c_str());
    string dQdyname = "dQdy" + to_string(i);
    dQdy[i] = op_decl_dat(mesh->cells, 15, "double", dQdy_data[i], dQdyname.c_str());
    string visBCname = "visBC" + to_string(i);
    visBC[i] = op_decl_dat(mesh->cells, 21, "double", visBC_data[i], visBCname.c_str());
  }
  divVelT   = op_decl_dat(mesh->cells, 15, "double", divVelT_data, "divVelT");
  curlVel   = op_decl_dat(mesh->cells, 15, "double", curlVel_data, "curlVel");
  pRHS      = op_decl_dat(mesh->cells, 15, "double", pRHS_data, "pRHS");
  pRHSex    = op_decl_dat(mesh->cells, 15, "double", pRHSex_data, "pRHSex");
  p         = op_decl_dat(mesh->cells, 15, "double", p_data, "p");
  dpdx      = op_decl_dat(mesh->cells, 15, "double", dpdx_data, "dpdx");
  dpdy      = op_decl_dat(mesh->cells, 15, "double", dpdy_data, "dpdy");
  prBC      = op_decl_dat(mesh->cells, 21, "double", prBC_data, "prBC");
  vorticity = op_decl_dat(mesh->cells, 15, "double", vorticity_data, "vorticity");
  save_temp = op_decl_dat(mesh->cells, 16, "double", save_temp_data, "save_temp");

  op_decl_const(1, "double", &gam);
  op_decl_const(1, "double", &mu);
  op_decl_const(1, "double", &nu);
  op_decl_const(1, "double", &bc_mach);
  op_decl_const(1, "double", &bc_alpha);
  op_decl_const(1, "double", &bc_p);
  op_decl_const(1, "double", &bc_u);
  op_decl_const(1, "double", &bc_v);
  op_decl_const(15, "int", FMASK);
  op_decl_const(1, "double", &ic_u);
  op_decl_const(1, "double", &ic_v);
  op_decl_const(46, "double", cubW_g);
  op_decl_const(46*15, "double", cubV_g);
  op_decl_const(46*15, "double", cubVDr_g);
  op_decl_const(46*15, "double", cubVDs_g);
  op_decl_const(7*15, "double", gF0Dr_g);
  op_decl_const(7*15, "double", gF0Ds_g);
  op_decl_const(7*15, "double", gF1Dr_g);
  op_decl_const(7*15, "double", gF1Ds_g);
  op_decl_const(7*15, "double", gF2Dr_g);
  op_decl_const(7*15, "double", gF2Ds_g);
  op_decl_const(7, "double", gaussW_g);
  op_decl_const(7*15, "double", gFInterp0_g);
  op_decl_const(7*15, "double", gFInterp1_g);
  op_decl_const(7*15, "double", gFInterp2_g);
  op_decl_const(7*15, "double", gF0DrR_g);
  op_decl_const(7*15, "double", gF0DsR_g);
  op_decl_const(7*15, "double", gF1DrR_g);
  op_decl_const(7*15, "double", gF1DsR_g);
  op_decl_const(7*15, "double", gF2DrR_g);
  op_decl_const(7*15, "double", gF2DsR_g);
  op_decl_const(7*15, "double", gFInterp0R_g);
  op_decl_const(7*15, "double", gFInterp1R_g);
  op_decl_const(7*15, "double", gFInterp2R_g);
}

INSData::~INSData() {
  for(int i = 0; i < 4; i++) {
    free(F_data[i]);
    free(div_data[i]);
  }
  for(int i = 0; i < 2; i++) {
    free(Q_data[0][i]);
    free(Q_data[1][i]);
    free(QT_data[i]);
    free(QTT_data[i]);

    free(N_data[0][i]);
    free(N_data[1][i]);
    free(exQ_data[i]);
    free(flux_data[i]);
    free(gradCurlVel_data[i]);
    free(dPdN_data[i]);
    free(visRHS_data[i]);
    free(visBC_data[i]);
    free(dQdx_data[i]);
    free(dQdy_data[i]);
  }
  free(divVelT_data);
  free(curlVel_data);
  free(pRHS_data);
  free(pRHSex_data);
  free(p_data);
  free(dpdx_data);
  free(dpdy_data);
  free(prBC_data);
  free(vorticity_data);
  free(save_temp_data);
}

void INSData::init() {

}

CubatureData::CubatureData(DGMesh *m, INSData *dat) {
  data = dat;
  mesh = m;

  Dx_data    = (double *)calloc(46 * 15 * mesh->numCells, sizeof(double));
  Dy_data    = (double *)calloc(46 * 15 * mesh->numCells, sizeof(double));
  OP_data    = (double *)calloc(15 * 15 * mesh->numCells, sizeof(double));
  temp_data  = (double *)calloc(46 * 15 * mesh->numCells, sizeof(double));
  temp2_data = (double *)calloc(46 * 15 * mesh->numCells, sizeof(double));

  Dx    = op_decl_dat(mesh->cells, 46 * 15, "double", Dx_data, "cub-Dx");
  Dy    = op_decl_dat(mesh->cells, 46 * 15, "double", Dy_data, "cub-Dy");
  OP    = op_decl_dat(mesh->cells, 15 * 15, "double", OP_data, "cub-OP");
  temp  = op_decl_dat(mesh->cells, 46 * 15, "double", temp_data, "cub-temp");
  temp2 = op_decl_dat(mesh->cells, 46 * 15, "double", temp2_data, "cub-temp2");
}

CubatureData::~CubatureData() {
  free(Dx_data);
  free(Dy_data);
  free(OP_data);
  free(temp_data);
  free(temp2_data);
}

void CubatureData::init() {
  // Initialise geometric factors for calcuating grad matrix
  op2_gemv(false, 46, 15, 1.0, constants->get_ptr(DGConstants::CUB_VDR), 15, mesh->x, 0.0, mesh->cubature->op_tmp[0]);
  op2_gemv(false, 46, 15, 1.0, constants->get_ptr(DGConstants::CUB_VDS), 15, mesh->x, 0.0, mesh->cubature->op_tmp[1]);
  op2_gemv(false, 46, 15, 1.0, constants->get_ptr(DGConstants::CUB_VDR), 15, mesh->y, 0.0, mesh->cubature->op_tmp[2]);
  op2_gemv(false, 46, 15, 1.0, constants->get_ptr(DGConstants::CUB_VDS), 15, mesh->y, 0.0, mesh->cubature->op_tmp[3]);

  op_par_loop(init_cubature_grad, "init_cubature_grad", mesh->cells,
              op_arg_dat(mesh->cubature->op_tmp[0], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[1], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[2], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(mesh->cubature->op_tmp[3], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(Dx, -1, OP_ID, 46 * 15, "double", OP_WRITE),
              op_arg_dat(Dy, -1, OP_ID, 46 * 15, "double", OP_WRITE));
  // Dx and Dy are row-major at this point

  // Calculate Cubature OP (contribution of Cubature points to Poisson matrix)
  op_par_loop(init_cubature_OP, "init_cubature_OP", mesh->cells,
              op_arg_dat(mesh->cubature->J, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(Dx,    -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(Dy,    -1, OP_ID, 46 * 15, "double", OP_READ),
              op_arg_dat(temp,  -1, OP_ID, 46 * 15, "double", OP_WRITE),
              op_arg_dat(temp2, -1, OP_ID, 46 * 15, "double", OP_WRITE));
  // Temp and temp2 are in row-major at this point
  op2_gemm_batch(false, true, 15, 15, 46, 1.0, Dx, 15, temp, 15, 0.0, OP, 15);
  op2_gemm_batch(false, true, 15, 15, 46, 1.0, Dy, 15, temp2, 15, 1.0, OP, 15);
  // OP is in col-major at this point
}

GaussData::GaussData(DGMesh *m, INSData *dat) {
  data = dat;
  mesh = m;

  rx_data      = (double *)calloc(21 * mesh->numCells, sizeof(double));
  sx_data      = (double *)calloc(21 * mesh->numCells, sizeof(double));
  ry_data      = (double *)calloc(21 * mesh->numCells, sizeof(double));
  sy_data      = (double *)calloc(21 * mesh->numCells, sizeof(double));
  tau_data     = (double *)calloc(3 * mesh->numCells, sizeof(double));
  reverse_data = (int *)calloc(3 * mesh->numCells, sizeof(int));
  for(int i = 0; i < 3; i++) {
    mDx_data[i] = (double *)calloc(7 * 15 * mesh->numCells, sizeof(double));
    mDy_data[i] = (double *)calloc(7 * 15 * mesh->numCells, sizeof(double));
    pDx_data[i] = (double *)calloc(7 * 15 * mesh->numCells, sizeof(double));
    pDy_data[i] = (double *)calloc(7 * 15 * mesh->numCells, sizeof(double));
    mD_data[i]  = (double *)calloc(7 * 15 * mesh->numCells, sizeof(double));
    pD_data[i]  = (double *)calloc(7 * 15 * mesh->numCells, sizeof(double));
    OP_data[i]  = (double *)calloc(15 * 15 * mesh->numCells, sizeof(double));
    OPf_data[i] = (double *)calloc(15 * 15 * mesh->numCells, sizeof(double));
  }

  rx      = op_decl_dat(mesh->cells, 21, "double", rx_data, "gauss-rx");
  sx      = op_decl_dat(mesh->cells, 21, "double", sx_data, "gauss-sx");
  ry      = op_decl_dat(mesh->cells, 21, "double", ry_data, "gauss-ry");
  sy      = op_decl_dat(mesh->cells, 21, "double", sy_data, "gauss-sy");
  tau     = op_decl_dat(mesh->cells, 3, "double", tau_data, "gauss-tau");
  reverse = op_decl_dat(mesh->cells, 3, "int", reverse_data, "gauss-reverse");
  for(int i = 0; i < 3; i++) {
    string name = "mDx" + to_string(i);
    mDx[i] = op_decl_dat(mesh->cells, 7 * 15, "double", mDx_data[i], name.c_str());
    name = "mDy" + to_string(i);
    mDy[i] = op_decl_dat(mesh->cells, 7 * 15, "double", mDy_data[i], name.c_str());
    name = "pDx" + to_string(i);
    pDx[i] = op_decl_dat(mesh->cells, 7 * 15, "double", pDx_data[i], name.c_str());
    name = "pDy" + to_string(i);
    pDy[i] = op_decl_dat(mesh->cells, 7 * 15, "double", pDy_data[i], name.c_str());
    name = "mD" + to_string(i);
    mD[i] = op_decl_dat(mesh->cells, 7 * 15, "double", mD_data[i], name.c_str());
    name = "pD" + to_string(i);
    pD[i] = op_decl_dat(mesh->cells, 7 * 15, "double", pD_data[i], name.c_str());
    name = "OP" + to_string(i);
    OP[i] = op_decl_dat(mesh->cells, 15 * 15, "double", OP_data[i], name.c_str());
    name = "OPf" + to_string(i);
    OPf[i] = op_decl_dat(mesh->cells, 15 * 15, "double", OPf_data[i], name.c_str());
  }
}

GaussData::~GaussData() {
  free(rx_data);
  free(sx_data);
  free(ry_data);
  free(sy_data);
  free(tau_data);
  free(reverse_data);
  for(int i = 0; i < 3; i++) {
    free(mDx_data[i]);
    free(mDy_data[i]);
    free(pDx_data[i]);
    free(pDy_data[i]);
    free(mD_data[i]);
    free(pD_data[i]);
    free(OP_data[i]);
    free(OPf_data[i]);
  }
}

void GaussData::init() {
  // Check which edges will require matrices to be 'reverse'
  op_par_loop(gauss_reverse, "gauss_reverse", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->nodeX,   -2, mesh->edge2cells, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY,   -2, mesh->edge2cells, 3, "double", OP_READ),
              op_arg_dat(reverse,       -2, mesh->edge2cells, 3, "int", OP_INC));

  // Calculate tau (used when constructing the Poisson matrix)
  op_par_loop(gauss_tau, "gauss_tau", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fscale,  -2, mesh->edge2cells, 15, "double", OP_READ),
              op_arg_dat(tau,           -2, mesh->edge2cells, 3, "double", OP_INC));

  op_par_loop(gauss_tau_bc, "gauss_tau_bc", mesh->bedges,
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->fscale,    0, mesh->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(tau,             0, mesh->bedge2cells, 3, "double", OP_INC));

  // Calculate geometric factors used when constructing gradient matrices
  init_gauss_grad_blas(mesh, this);

  op_par_loop(init_gauss_grad, "init_gauss_grad", mesh->cells,
              op_arg_dat(rx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(ry, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sy, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(mDx[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[2], -1, OP_ID, 7 * 15, "double", OP_WRITE));

  // Construct gradient matrices (1 for each face)
  op_par_loop(init_gauss_grad2, "init_gauss_grad2", mesh->cells,
              op_arg_dat(mesh->gauss->nx, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(mDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mD[0],  -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mD[1],  -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mD[2],  -1, OP_ID, 7 * 15, "double", OP_WRITE));

  // Calculate geometric factors for grad matrices used by neighbours
  // Matrices are calculated locally, then copied to neighbour elements
  init_gauss_grad_neighbour_blas(mesh, this);

  op_par_loop(init_gauss_grad_neighbour, "init_gauss_grad_neighbour", mesh->cells,
              op_arg_dat(reverse, -1, OP_ID, 3, "int", OP_READ),
              op_arg_dat(rx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sx, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(ry, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(sy, -1, OP_ID, 21, "double", OP_RW),
              op_arg_dat(mDx[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[2], -1, OP_ID, 7 * 15, "double", OP_WRITE));

  // Copy x and y grad matrices to neighbours
  op_par_loop(gauss_grad_faces, "gauss_grad_faces", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mDx[0], -2, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[0], -2, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[1], -2, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[1], -2, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[2], -2, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[2], -2, mesh->edge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(pDx[0], -2, mesh->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[0], -2, mesh->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDx[1], -2, mesh->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[1], -2, mesh->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDx[2], -2, mesh->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[2], -2, mesh->edge2cells, 7 * 15, "double", OP_INC));

  // Calculate final neighbour grad matrix
  op_par_loop(init_gauss_grad2, "init_gauss_grad2", mesh->cells,
              op_arg_dat(mesh->gauss->nx, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -1, OP_ID, 21, "double", OP_READ),
              op_arg_dat(pDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDy[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDy[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDx[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDy[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pD[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pD[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pD[2], -1, OP_ID, 7 * 15, "double", OP_WRITE));

  // Calculate Gauss OP for each face (local contribution of face in Poisson matrix)
  // Face 0 temps: mDx, Face 1 temps: mDy, Face 2 temps: pDx
  op_par_loop(gauss_op, "gauss_op", mesh->cells,
              op_arg_dat(tau, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -1, OP_ID, 21, "double", OP_READ),
              // Face 0
              op_arg_dat(mD[0], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDx[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDx[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              // Face 1
              op_arg_dat(mD[1], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(mDy[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(mDy[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              // Face 2
              op_arg_dat(mD[2], -1, OP_ID, 7 * 15, "double", OP_READ),
              op_arg_dat(pDx[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pDx[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pDx[2], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              // Reset dats for OPf
              op_arg_dat(pDy[0], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pDy[1], -1, OP_ID, 7 * 15, "double", OP_WRITE),
              op_arg_dat(pDy[2], -1, OP_ID, 7 * 15, "double", OP_WRITE));

  op2_gemm(true, false, 15, 15, 7, 1.0, mDx[0], 7, constants->get_ptr(DGConstants::GAUSS_FINTERP0), 7, 0.0, OP[0], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, mDx[1], 7, mD[0], 15, 1.0, OP[0], 15);
  op2_gemm(true, false, 15, 15, 7, -1.0, mDx[2], 7, constants->get_ptr(DGConstants::GAUSS_FINTERP0), 7, 1.0, OP[0], 15);

  op2_gemm(true, false, 15, 15, 7, 1.0, mDy[0], 7, constants->get_ptr(DGConstants::GAUSS_FINTERP1), 7, 0.0, OP[1], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, mDy[1], 7, mD[1], 15, 1.0, OP[1], 15);
  op2_gemm(true, false, 15, 15, 7, -1.0, mDy[2], 7, constants->get_ptr(DGConstants::GAUSS_FINTERP1), 7, 1.0, OP[1], 15);

  op2_gemm(true, false, 15, 15, 7, 1.0, pDx[0], 7, constants->get_ptr(DGConstants::GAUSS_FINTERP2), 7, 0.0, OP[2], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, pDx[1], 7, mD[2], 15, 1.0, OP[2], 15);
  op2_gemm(true, false, 15, 15, 7, -1.0, pDx[2], 7, constants->get_ptr(DGConstants::GAUSS_FINTERP2), 7, 1.0, OP[2], 15);

  // Calculate Gauss OPf for each face (contribution to neighbouring element in Poisson matrix)
  op_par_loop(gauss_gfi_faces, "gauss_gfi_faces", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(pDy[0], -2, mesh->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[1], -2, mesh->edge2cells, 7 * 15, "double", OP_INC),
              op_arg_dat(pDy[2], -2, mesh->edge2cells, 7 * 15, "double", OP_INC));

  op2_gemm_batch(true, true, 15, 15, 7, 1.0, mDx[0], 7, pDy[0], 15, 0.0, OPf[0], 15);
  op2_gemm_batch(true, true, 15, 15, 7, 1.0, mDx[1], 7, pD[0], 15, 1.0, OPf[0], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, mDx[2], 7, pDy[0], 15, 1.0, OPf[0], 15);

  op2_gemm_batch(true, true, 15, 15, 7, 1.0, mDy[0], 7, pDy[1], 15, 0.0, OPf[1], 15);
  op2_gemm_batch(true, true, 15, 15, 7, 1.0, mDy[1], 7, pD[1], 15, 1.0, OPf[1], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, mDy[2], 7, pDy[1], 15, 1.0, OPf[1], 15);

  op2_gemm_batch(true, true, 15, 15, 7, 1.0, pDx[0], 7, pDy[2], 15, 0.0, OPf[2], 15);
  op2_gemm_batch(true, true, 15, 15, 7, 1.0, pDx[1], 7, pD[2], 15, 1.0, OPf[2], 15);
  op2_gemm_batch(true, true, 15, 15, 7, -1.0, pDx[2], 7, pDy[2], 15, 1.0, OPf[2], 15);

  // Applying the correct factors to OP and OPf is done when constructing the Poisson matrix
}
