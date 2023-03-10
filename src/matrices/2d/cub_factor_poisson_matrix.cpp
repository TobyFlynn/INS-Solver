#include "matrices/2d/cub_factor_poisson_matrix_2d.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"

extern Timing *timer;
extern DGConstants *constants;

CubFactorPoissonMatrix2D::CubFactorPoissonMatrix2D(DGMesh2D *m) : FactorPoissonMatrix2D(m) {
  DG_FP *tmp_cub_np = (DG_FP *)calloc(DG_CUB_NP * mesh->cells->size, sizeof(DG_FP));
  cFactor  = op_decl_dat(mesh->cells, DG_CUB_NP, DG_FP_STR, tmp_cub_np, "poisson_cFactor");
  free(tmp_cub_np);
}

void CubFactorPoissonMatrix2D::calc_op1() {
  timer->startTimer("CubFactorPoissonMatrix2D - calc_op1");
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
  timer->endTimer("CubFactorPoissonMatrix2D - calc_op1");
}
