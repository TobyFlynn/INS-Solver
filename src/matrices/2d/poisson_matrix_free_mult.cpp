#include "matrices/2d/poisson_matrix_free_mult_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"
#include "dg_dat_pool.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;
extern DGDatPool *dg_dat_pool;

PoissonMatrixFreeMult2D::PoissonMatrixFreeMult2D(DGMesh2D *m) {
  mesh = m;
}

void PoissonMatrixFreeMult2D::mat_free_set_bc_types(op_dat bc_ty) {
  mat_free_bcs = bc_ty;
}

void PoissonMatrixFreeMult2D::mat_free_apply_bc(op_dat rhs, op_dat bc) {
  if(mesh->bface2cells) {
    // TODO
  }
}

void PoissonMatrixFreeMult2D::mat_free_mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonMatrixFreeMult2D - mult");
  DGTempDat tmp_grad0 = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_grad1 = dg_dat_pool->requestTempDatCells(DG_NP);
  timer->startTimer("PoissonMatrixFreeMult2D - mult grad");
  mesh->grad(in, tmp_grad0.dat, tmp_grad1.dat);
  timer->endTimer("PoissonMatrixFreeMult2D - mult grad");

  DGTempDat tmp_npf0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_npf1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_npf2 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(tmp_npf0.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf1.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf2.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("PoissonMatrixFreeMult2D - mult faces");
  op_par_loop(pmf_2d_mult_faces, "pmf_2d_mult_faces", mesh->faces,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale,  -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ,  -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(in, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_grad0.dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_grad1.dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf0.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf1.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf2.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonMatrixFreeMult2D - mult faces");

  timer->startTimer("PoissonMatrixFreeMult2D - mult bfaces");
  if(mesh->bface2cells) {
    // TODO
  }
  timer->endTimer("PoissonMatrixFreeMult2D - mult bfaces");

  timer->startTimer("PoissonMatrixFreeMult2D - mult cells");
  timer->startTimer("PoissonMatrixFreeMult2D - mult cells MM");
  mesh->mass(tmp_grad0.dat);
  mesh->mass(tmp_grad1.dat);
  timer->endTimer("PoissonMatrixFreeMult2D - mult cells MM");

  timer->startTimer("PoissonMatrixFreeMult2D - mult cells Emat");
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf0.dat, 1.0, tmp_grad0.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf1.dat, 1.0, tmp_grad1.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf2.dat, 0.0, out);
  timer->endTimer("PoissonMatrixFreeMult2D - mult cells Emat");

  timer->startTimer("PoissonMatrixFreeMult2D - mult cells cells");
  op_par_loop(pmf_2d_mult_cells_geof, "pmf_2d_mult_cells_geof", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_grad0.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_grad1.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, true, 1.0, DGConstants::DR, tmp_grad0.dat, 1.0, out);
  op2_gemv(mesh, true, 1.0, DGConstants::DS, tmp_grad1.dat, 1.0, out);
  timer->endTimer("PoissonMatrixFreeMult2D - mult cells cells");
  timer->endTimer("PoissonMatrixFreeMult2D - mult cells");
  dg_dat_pool->releaseTempDatCells(tmp_grad0);
  dg_dat_pool->releaseTempDatCells(tmp_grad1);
  dg_dat_pool->releaseTempDatCells(tmp_npf0);
  dg_dat_pool->releaseTempDatCells(tmp_npf1);
  dg_dat_pool->releaseTempDatCells(tmp_npf2);
  timer->endTimer("PoissonMatrixFreeMult2D - mult");
}
