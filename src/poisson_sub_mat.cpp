#include "poisson.h"

#include "op_seq.h"

#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_global_constants.h"

void PoissonSolve::calc_cub_sub_mat() {
  // Initialise geometric factors for calcuating grad matrix
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDR, mesh->x, 0.0, mesh->cubature->op_tmp[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDS, mesh->x, 0.0, mesh->cubature->op_tmp[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDR, mesh->y, 0.0, mesh->cubature->op_tmp[2]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDS, mesh->y, 0.0, mesh->cubature->op_tmp[3]);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB_V, factor, 0.0, cFactor);

  op_par_loop(poisson_cubature_op, "poisson_cubature_op", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(cubVDr_g, DG_ORDER * DG_CUB_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(cubVDs_g, DG_ORDER * DG_CUB_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[0], -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[1], -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[2], -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[3], -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->J, -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(cFactor, -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_WRITE));
}

void PoissonSolve::calc_gauss_sub_mat() {
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, factor, 0.0, gFactor);

  // TODO change for edges with different orders on each side
  op_par_loop(poisson_gauss_grad, "poisson_gauss_grad", mesh->edges,
              op_arg_dat(mesh->order, -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_gbl(gF0Dr_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF0Ds_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF1Dr_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF1Ds_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF2Dr_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF2Ds_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gFInterp0_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gFInterp1_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gFInterp2_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->x, -2, mesh->edge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -2, mesh->edge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(h,       -2, mesh->edge2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->fscale, -2, mesh->edge2cells, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(op1, 0, mesh->edge2cells, DG_NP * DG_NP, "double", OP_INC),
              op_arg_dat(op1, 1, mesh->edge2cells, DG_NP * DG_NP, "double", OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_WRITE));

  // If not dirichlet BC, kernel will assume it is a neumann bc
  op_par_loop(poisson_gauss_grad_b, "poisson_gauss_grad_b", mesh->bedges,
              op_arg_dat(mesh->order, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_gbl(gF0Dr_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF0Ds_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF1Dr_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF1Ds_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF2Dr_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gF2Ds_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gFInterp0_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gFInterp1_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(gFInterp2_g, DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(h,       0, mesh->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->fscale, 0, mesh->bedge2cells, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(op1, 0, mesh->bedge2cells, DG_NP * DG_NP, "double", OP_INC),
              op_arg_dat(op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE));
}

void PoissonSolve::calc_mm_mat() {
  op_par_loop(poisson_mm, "poisson_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mmFactor, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_INC));
}
