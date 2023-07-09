#include "matrices/2d/factor_poisson_matrix_free_mult_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"
#include "dg_dat_pool.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;
extern DGDatPool *dg_dat_pool;

FactorPoissonMatrixFreeMult2D::FactorPoissonMatrixFreeMult2D(DGMesh2D *m) : PoissonMatrixFreeMult2D(m) {
  mat_free_gtau = op_decl_dat(mesh->cells, 3, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_free_tau");
  mat_free_factor_copy = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_free_factor_copy");
}

void FactorPoissonMatrixFreeMult2D::calc_tau() {
  timer->startTimer("FactorPoissonMatrixFreeMult2D - calc tau");
  op_par_loop(fpmf_2d_calc_tau_faces, "fpmf_2d_calc_tau_faces", mesh->faces,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor_copy, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_gtau, -2, mesh->face2cells, 3, DG_FP_STR, OP_WRITE));
  if(mesh->bface2cells) {
    op_par_loop(fpmf_2d_calc_tau_bfaces, "fpmf_2d_calc_tau_bfaces", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_factor_copy, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_gtau, 0, mesh->bface2cells, 3, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("FactorPoissonMatrixFreeMult2D - calc tau");
}

void FactorPoissonMatrixFreeMult2D::mat_free_set_factor(op_dat f) {
  mat_free_factor = f;

  factor_order = mesh->order_int;
  current_order = mesh->order_int;

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(mat_free_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor_copy, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  calc_tau();
}

void FactorPoissonMatrixFreeMult2D::check_current_order() {
  timer->startTimer("FactorPoissonMatrixFreeMult2D - check order");
  if(current_order != mesh->order_int) {
    mesh->interp_dat_between_orders(current_order, mesh->order_int, mat_free_factor, mat_free_factor_copy);
    current_order = mesh->order_int;

    calc_tau();
  }
  timer->endTimer("FactorPoissonMatrixFreeMult2D - check order");
}

void FactorPoissonMatrixFreeMult2D::mat_free_apply_bc(op_dat rhs, op_dat bc) {
  if(mesh->bface2cells) {
    check_current_order();
    op_par_loop(fpmf_2d_apply_bc, "fpmf_2d_apply_bc", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mat_free_bcs, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, 0, mesh->bface2cells, 5, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_factor_copy, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(rhs, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
}

void FactorPoissonMatrixFreeMult2D::mat_free_mult(op_dat in, op_dat out) {
  check_current_order();
  timer->startTimer("FactorPoissonMatrixFreeMult2D - mult");
  DGTempDat tmp_grad0 = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_grad1 = dg_dat_pool->requestTempDatCells(DG_NP);
  timer->startTimer("FactorPoissonMatrixFreeMult2D - mult grad");
  op2_gemv(mesh, false, 1.0, DGConstants::DR, in, 0.0, tmp_grad0.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::DS, in, 0.0, tmp_grad1.dat);
  op_par_loop(fpmf_2d_grad, "fpmf_2d_grad", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor_copy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_grad0.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_grad1.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("FactorPoissonMatrixFreeMult2D - mult grad");

  DGTempDat tmp_npf0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_npf1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_npf2 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(tmp_npf0.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf1.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf2.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("FactorPoissonMatrixFreeMult2D - mult faces");
  mesh->jump(in, tmp_npf0.dat);
  mesh->avg(tmp_grad0.dat, tmp_npf1.dat);
  mesh->avg(tmp_grad1.dat, tmp_npf2.dat);
  timer->endTimer("FactorPoissonMatrixFreeMult2D - mult faces");

  timer->startTimer("FactorPoissonMatrixFreeMult2D - mult bfaces");
  if(mesh->bface2cells) {
    op_par_loop(pmf_2d_mult_avg_jump, "pmf_2d_mult_avg_jump", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mat_free_bcs, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_grad0.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_grad1.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf0.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf1.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf2.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }
  timer->endTimer("FactorPoissonMatrixFreeMult2D - mult bfaces");

  timer->startTimer("FactorPoissonMatrixFreeMult2D - finish flux");
    op_par_loop(fpmf_2d_mult_flux, "fpmf_2d_mult_flux", mesh->cells,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mesh->nx_c_new, -1, OP_ID, 3, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->ny_c_new, -1, OP_ID, 3, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sJ_c_new, -1, OP_ID, 3, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_gtau, -1, OP_ID, 3, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_factor_copy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf0.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(tmp_npf1.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(tmp_npf2.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));
  timer->endTimer("FactorPoissonMatrixFreeMult2D - finish flux");

  timer->startTimer("FactorPoissonMatrixFreeMult2D - mult cells");
  timer->startTimer("FactorPoissonMatrixFreeMult2D - mult cells MM");
  mesh->mass(tmp_grad0.dat);
  mesh->mass(tmp_grad1.dat);
  timer->endTimer("FactorPoissonMatrixFreeMult2D - mult cells MM");

  timer->startTimer("FactorPoissonMatrixFreeMult2D - mult cells Emat");
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf1.dat, 1.0, tmp_grad0.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf2.dat, 1.0, tmp_grad1.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf0.dat, 0.0, out);
  timer->endTimer("FactorPoissonMatrixFreeMult2D - mult cells Emat");

  timer->startTimer("FactorPoissonMatrixFreeMult2D - mult cells cells");
  op_par_loop(pmf_2d_mult_cells_geof, "pmf_2d_mult_cells_geof", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_grad0.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_grad1.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, true, 1.0, DGConstants::DR, tmp_grad0.dat, 1.0, out);
  op2_gemv(mesh, true, 1.0, DGConstants::DS, tmp_grad1.dat, 1.0, out);
  timer->endTimer("FactorPoissonMatrixFreeMult2D - mult cells cells");
  timer->endTimer("FactorPoissonMatrixFreeMult2D - mult cells");
  dg_dat_pool->releaseTempDatCells(tmp_grad0);
  dg_dat_pool->releaseTempDatCells(tmp_grad1);
  dg_dat_pool->releaseTempDatCells(tmp_npf0);
  dg_dat_pool->releaseTempDatCells(tmp_npf1);
  dg_dat_pool->releaseTempDatCells(tmp_npf2);
  timer->endTimer("FactorPoissonMatrixFreeMult2D - mult");
}
