#include "matrices/2d/factor_poisson_matrix_2d.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"

extern Timing *timer;
extern DGConstants *constants;

FactorPoissonMatrix2D::FactorPoissonMatrix2D(DGMesh2D *m) : PoissonMatrix2D(m) {
  DG_FP *tmp_g_np = (DG_FP *)calloc(DG_G_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *tmp_cub_np = (DG_FP *)calloc(DG_CUB_NP * mesh->cells->size, sizeof(DG_FP));
  gFactor  = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np, "poisson_gFactor");
  gS       = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np, "poisson_gS");
  cFactor  = op_decl_dat(mesh->cells, DG_CUB_NP, DG_FP_STR, tmp_cub_np, "poisson_cFactor");
  free(tmp_cub_np);
  free(tmp_g_np);
}

void FactorPoissonMatrix2D::set_factor(op_dat f) {
  factor = f;
}

void FactorPoissonMatrix2D::set_surface(op_dat surf, const DG_FP al) {
  s = surf;
  alpha = al;
}

void FactorPoissonMatrix2D::calc_op1() {
  timer->startTimer("FactorPoissonMatrix2D - calc_op1");
  // Initialise geometric factors for calcuating grad matrix
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDR, mesh->x, 0.0, mesh->cubature->op_tmp[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDS, mesh->x, 0.0, mesh->cubature->op_tmp[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDR, mesh->y, 0.0, mesh->cubature->op_tmp[2]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDS, mesh->y, 0.0, mesh->cubature->op_tmp[3]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_V, factor, 0.0, cFactor);

  op_par_loop(fact_poisson_cub_op1, "fact_poisson_cub_op1", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::CUB_VDR), DG_ORDER * DG_CUB_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::CUB_VDS), DG_ORDER * DG_CUB_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[0], -1, OP_ID, DG_CUB_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[1], -1, OP_ID, DG_CUB_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[2], -1, OP_ID, DG_CUB_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[3], -1, OP_ID, DG_CUB_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->cubature->J, -1, OP_ID, DG_CUB_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cFactor, -1, OP_ID, DG_CUB_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("FactorPoissonMatrix2D - calc_op1");
}

void FactorPoissonMatrix2D::calc_op2() {
  timer->startTimer("FactorPoissonMatrix2D - calc_op2");
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, factor, 0.0, gFactor);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, s, 0.0, gS);

  op_par_loop(fact_poisson_gauss_op2, "fact_poisson_gauss_op2", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP0), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP1), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP2), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ), 
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->x, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -2, mesh->face2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(gFactor, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gS, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, 0, mesh->face2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC),
              op_arg_dat(op1, 1, mesh->face2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("FactorPoissonMatrix2D - calc_op2");
}

void FactorPoissonMatrix2D::calc_opbc() {
  timer->startTimer("FactorPoissonMatrix2D - calc_opbc");
  if(mesh->bface2cells) {
    op_par_loop(fact_poisson_gauss_bop, "fact_poisson_gauss_bop", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP0), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP1), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP2), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->fscale, 0, mesh->bface2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(gFactor, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(op1, 0, mesh->bface2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC),
                op_arg_dat(opbc, -1, OP_ID, DG_GF_NP * DG_NP, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("FactorPoissonMatrix2D - calc_opbc");
}
