#include "matrices/2d/factor_poisson_coarse_matrix_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

FactorPoissonCoarseMatrix2D::FactorPoissonCoarseMatrix2D(DGMesh2D *m) : PoissonCoarseMatrix2D(m) {
  factor = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_coarse_factor");
}

void FactorPoissonCoarseMatrix2D::set_factor(op_dat f) {
  mesh->interp_dat_between_orders(DG_ORDER, 1, f, factor);
}

void FactorPoissonCoarseMatrix2D::calc_op1() {
  timer->startTimer("FactorPoissonCoarseMatrix2D - calc_op1");
  op_par_loop(factor_poisson_coarse_matrix_2d_op1, "factor_poisson_coarse_matrix_2d_op1", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(factor,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE));
  timer->endTimer("FactorPoissonCoarseMatrix2D - calc_op1");
}

void FactorPoissonCoarseMatrix2D::calc_op2() {
  timer->startTimer("FactorPoissonCoarseMatrix2D - calc_op2");
  op_par_loop(factor_poisson_coarse_matrix_2d_op2, "factor_poisson_coarse_matrix_2d_op2", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -2, mesh->face2cells, 5, DG_FP_STR, OP_READ),
              op_arg_dat(factor, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, 0, mesh->face2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC),
              op_arg_dat(op1, 1, mesh->face2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE));
  timer->endTimer("FactorPoissonCoarseMatrix2D - calc_op2");
}

void FactorPoissonCoarseMatrix2D::calc_opbc() {
  timer->startTimer("FactorPoissonCoarseMatrix2D - calc_opbc");
  if(mesh->bface2cells) {
    op_par_loop(factor_poisson_coarse_matrix_2d_bop, "factor_poisson_coarse_matrix_2d_bop", mesh->bfaces,
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, 0, mesh->bface2cells, 5, DG_FP_STR, OP_READ),
                op_arg_dat(factor,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(op1, 0, mesh->bface2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC));

    op_par_loop(factor_poisson_coarse_matrix_2d_opbc, "factor_poisson_coarse_matrix_2d_opbc", mesh->bfaces,
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, 0, mesh->bface2cells, 5, DG_FP_STR, OP_READ),
                op_arg_dat(factor,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(opbc, -1, OP_ID, DG_NPF_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("FactorPoissonCoarseMatrix2D - calc_opbc");
}
