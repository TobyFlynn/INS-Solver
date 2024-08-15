#include "solvers/2d/compressible_euler.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "dg_utils.h"
#include "dg_abort.h"
#include "timing.h"
#include "config.h"

#include <string>

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

CompressibleEuler2D::CompressibleEuler2D(DGMesh2D *m) {
  mesh = m;

  std::string name;
  for(int i = 0; i < 4; i++) {
    name = "Q" + std::to_string(i);
    Q[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "rk_wQ" + std::to_string(i);
    rk_wQ[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "rk_RHSQ0" + std::to_string(i);
    rk_RHSQ[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "rk_RHSQ1" + std::to_string(i);
    rk_RHSQ[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "rk_RHSQ2" + std::to_string(i);
    rk_RHSQ[2][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
  }

  time = 0.0;
  over_int = true;
}

CompressibleEuler2D::~CompressibleEuler2D() {

}

void CompressibleEuler2D::init() {
  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;

  op_par_loop(euler_2d_ic, "euler_2d_ic", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, "double", OP_WRITE));
}

void CompressibleEuler2D::step() {
  if(over_int)
    rhs_over_int(Q, rk_RHSQ[0]);
  else
    rhs(Q, rk_RHSQ[0]);

  op_par_loop(euler_2d_wQ_0, "euler_2d_wQ_0", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_wQ[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  if(over_int)
    rhs_over_int(rk_wQ, rk_RHSQ[1]);
  else
    rhs(rk_wQ, rk_RHSQ[1]);

  op_par_loop(euler_2d_wQ_1, "euler_2d_wQ_1", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_wQ[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  if(over_int)
    rhs_over_int(rk_wQ, rk_RHSQ[2]);
  else
    rhs(rk_wQ, rk_RHSQ[2]);

  op_par_loop(euler_2d_wQ_2, "euler_2d_wQ_2", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[2][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[2][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[2][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[2][3], -1, OP_ID, DG_NP, "double", OP_READ));

  time += dt;
}

void CompressibleEuler2D::rhs(op_dat *inQ, op_dat *outQ) {
  DGTempDat F[4], G[4];
  F[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  F[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  F[2] = dg_dat_pool->requestTempDatCells(DG_NP);
  F[3] = dg_dat_pool->requestTempDatCells(DG_NP);
  G[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  G[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  G[2] = dg_dat_pool->requestTempDatCells(DG_NP);
  G[3] = dg_dat_pool->requestTempDatCells(DG_NP);

  timer->startTimer("Euler2D - Vol Kernel");
  op_par_loop(euler_2d_vol, "euler_2d_vol", mesh->cells,
              op_arg_dat(inQ[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[3], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(F[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(F[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(F[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(F[3].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[3].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("Euler2D - Vol Kernel");

  timer->startTimer("Euler2D - Vol Div");
  mesh->div_weak(F[0].dat, G[0].dat, outQ[0]);
  mesh->div_weak(F[1].dat, G[1].dat, outQ[1]);
  mesh->div_weak(F[2].dat, G[2].dat, outQ[2]);
  mesh->div_weak(F[3].dat, G[3].dat, outQ[3]);
  timer->endTimer("Euler2D - Vol Div");

  dg_dat_pool->releaseTempDatCells(F[0]);
  dg_dat_pool->releaseTempDatCells(F[1]);
  dg_dat_pool->releaseTempDatCells(F[2]);
  dg_dat_pool->releaseTempDatCells(F[3]);
  dg_dat_pool->releaseTempDatCells(G[0]);
  dg_dat_pool->releaseTempDatCells(G[1]);
  dg_dat_pool->releaseTempDatCells(G[2]);
  dg_dat_pool->releaseTempDatCells(G[3]);

  DGTempDat surf_terms[4];
  surf_terms[0] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  surf_terms[1] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  surf_terms[2] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  surf_terms[3] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(surf_terms[0].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(surf_terms[1].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(surf_terms[2].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(surf_terms[3].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("Euler2D - Surf Kernel");
  op_par_loop(euler_2d_surf, "euler_2d_surf", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[2], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[3], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(surf_terms[0].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(surf_terms[1].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(surf_terms[2].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(surf_terms[3].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  timer->endTimer("Euler2D - Surf Kernel");

  timer->startTimer("Euler2D - BC Surf Kernel");
  if(mesh->bface2cells) {
    op_par_loop(euler_2d_surf_bc, "euler_2d_surf_bc", mesh->bfaces,
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(inQ[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(inQ[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(inQ[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(inQ[3], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(surf_terms[0].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(surf_terms[1].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(surf_terms[2].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(surf_terms[3].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }
  timer->endTimer("Euler2D - BC Surf Kernel");

  timer->startTimer("Euler2D - Surf Lift");
  op2_gemv(mesh, false, -1.0, DGConstants::LIFT, surf_terms[0].dat, 1.0, outQ[0]);
  op2_gemv(mesh, false, -1.0, DGConstants::LIFT, surf_terms[1].dat, 1.0, outQ[1]);
  op2_gemv(mesh, false, -1.0, DGConstants::LIFT, surf_terms[2].dat, 1.0, outQ[2]);
  op2_gemv(mesh, false, -1.0, DGConstants::LIFT, surf_terms[3].dat, 1.0, outQ[3]);
  timer->endTimer("Euler2D - Surf Lift");

  dg_dat_pool->releaseTempDatCells(surf_terms[0]);
  dg_dat_pool->releaseTempDatCells(surf_terms[1]);
  dg_dat_pool->releaseTempDatCells(surf_terms[2]);
  dg_dat_pool->releaseTempDatCells(surf_terms[3]);
}

void CompressibleEuler2D::rhs_over_int(op_dat *inQ, op_dat *outQ) {
  DGTempDat F[4], G[4];
  F[0] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  F[1] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  F[2] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  F[3] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  G[0] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  G[1] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  G[2] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  G[3] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);

  timer->startTimer("Euler2D - Vol Interp");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, inQ[0], 0.0, F[0].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, inQ[1], 0.0, F[1].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, inQ[2], 0.0, F[2].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, inQ[3], 0.0, F[3].dat);
  timer->endTimer("Euler2D - Vol Interp");

  timer->startTimer("Euler2D - Vol Kernel");
  op_par_loop(euler_2d_vol_oi, "euler_2d_vol_oi", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(F[0].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(F[1].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(F[2].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(F[3].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(G[0].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[1].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[2].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(G[3].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("Euler2D - Vol Kernel");

  timer->startTimer("Euler2D - Vol Div");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, F[0].dat, 0.0, outQ[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, G[0].dat, 1.0, outQ[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, F[1].dat, 0.0, outQ[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, G[1].dat, 1.0, outQ[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, F[2].dat, 0.0, outQ[2]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, G[2].dat, 1.0, outQ[2]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, F[3].dat, 0.0, outQ[3]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, G[3].dat, 1.0, outQ[3]);
  timer->endTimer("Euler2D - Vol Div");

  dg_dat_pool->releaseTempDatCells(F[0]);
  dg_dat_pool->releaseTempDatCells(F[1]);
  dg_dat_pool->releaseTempDatCells(F[2]);
  dg_dat_pool->releaseTempDatCells(F[3]);
  dg_dat_pool->releaseTempDatCells(G[0]);
  dg_dat_pool->releaseTempDatCells(G[1]);
  dg_dat_pool->releaseTempDatCells(G[2]);
  dg_dat_pool->releaseTempDatCells(G[3]);

  DGTempDat surf_terms[4];
  surf_terms[0] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  surf_terms[1] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  surf_terms[2] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  surf_terms[3] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);

  op_par_loop(zero_cub_surf_2d, "zero_cub_surf_2d", mesh->cells,
                op_arg_dat(surf_terms[0].dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_cub_surf_2d, "zero_cub_surf_2d", mesh->cells,
                op_arg_dat(surf_terms[1].dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_cub_surf_2d, "zero_cub_surf_2d", mesh->cells,
                op_arg_dat(surf_terms[2].dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_cub_surf_2d, "zero_cub_surf_2d", mesh->cells,
                op_arg_dat(surf_terms[3].dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE));

  timer->startTimer("Euler2D - Surf Kernel");
  op_par_loop(euler_2d_surf_oi, "euler_2d_surf_oi", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[2], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(inQ[3], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(surf_terms[0].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(surf_terms[1].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(surf_terms[2].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(surf_terms[3].dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("Euler2D - Surf Kernel");

  timer->startTimer("Euler2D - BC Surf Kernel");
  if(mesh->bface2cells) {
    op_par_loop(euler_2d_surf_bc_oi, "euler_2d_surf_bc_oi", mesh->bfaces,
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(inQ[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(inQ[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(inQ[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(inQ[3], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(surf_terms[0].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC),
                op_arg_dat(surf_terms[1].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC),
                op_arg_dat(surf_terms[2].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC),
                op_arg_dat(surf_terms[3].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC));
  }
  timer->endTimer("Euler2D - BC Surf Kernel");

  timer->startTimer("Euler2D - Surf Lift");
  op2_gemv(mesh, false, -1.0, DGConstants::CUBSURF2D_LIFT, surf_terms[0].dat, 1.0, outQ[0]);
  op2_gemv(mesh, false, -1.0, DGConstants::CUBSURF2D_LIFT, surf_terms[1].dat, 1.0, outQ[1]);
  op2_gemv(mesh, false, -1.0, DGConstants::CUBSURF2D_LIFT, surf_terms[2].dat, 1.0, outQ[2]);
  op2_gemv(mesh, false, -1.0, DGConstants::CUBSURF2D_LIFT, surf_terms[3].dat, 1.0, outQ[3]);
  timer->endTimer("Euler2D - Surf Lift");

  dg_dat_pool->releaseTempDatCells(surf_terms[0]);
  dg_dat_pool->releaseTempDatCells(surf_terms[1]);
  dg_dat_pool->releaseTempDatCells(surf_terms[2]);
  dg_dat_pool->releaseTempDatCells(surf_terms[3]);
}

DG_FP CompressibleEuler2D::get_time() {
  return time;
}

void CompressibleEuler2D::set_dt(const DG_FP _dt) {
  dt = _dt;
}