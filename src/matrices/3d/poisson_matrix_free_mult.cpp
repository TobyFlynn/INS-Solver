#include "matrices/3d/poisson_matrix_free_mult_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"
#include "dg_dat_pool.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;
extern DGDatPool3D *dg_dat_pool;

PoissonMatrixFreeMult3D::PoissonMatrixFreeMult3D(DGMesh3D *m) {
  mesh = m;
  mat_free_tau_c = op_decl_dat(mesh->cells, 4, DG_FP_STR, (DG_FP *)NULL, "mat_free_tau_c");
}

void PoissonMatrixFreeMult3D::mat_free_set_bc_types(op_dat bc_ty) {
  mat_free_bcs = bc_ty;
}

void PoissonMatrixFreeMult3D::mat_free_apply_bc(op_dat rhs, op_dat bc) {
  if(mesh->bface2cells) {
    op_par_loop(pmf_3d_apply_bc, "pmf_3d_apply_bc", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mat_free_bcs, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, 0, mesh->bface2cells, 10, DG_FP_STR, OP_READ),
                op_arg_dat(bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(rhs, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
}

void PoissonMatrixFreeMult3D::mat_free_mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonMatrixFreeMult3D - mult");
  timer->startTimer("PoissonMatrixFreeMult3D - calc tau");
  op_par_loop(pmf_3d_calc_tau_faces, "pmf_3d_calc_tau_faces", mesh->faces,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_tau_c, -2, mesh->face2cells, 4, DG_FP_STR, OP_WRITE));
  if(mesh->bface2cells) {
    op_par_loop(pmf_3d_calc_tau_bfaces, "pmf_3d_calc_tau_bfaces", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_tau_c, 0, mesh->bface2cells, 4, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("PoissonMatrixFreeMult3D - calc tau");

  DGTempDat tmp_grad0 = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_grad1 = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_grad2 = dg_dat_pool->requestTempDatCells(DG_NP);
  timer->startTimer("PoissonMatrixFreeMult3D - mult grad");
  // mesh->grad(in, tmp_grad0.dat, tmp_grad1.dat, tmp_grad2.dat);
  mesh->grad_halo_exchange(in, tmp_grad0.dat, tmp_grad1.dat, tmp_grad2.dat);
  timer->endTimer("PoissonMatrixFreeMult3D - mult grad");

  DGTempDat tmp_npf0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_npf1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_npf2 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_npf3 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(tmp_npf0.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(tmp_npf1.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf2.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf3.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("PoissonMatrixFreeMult3D - mult faces");
  mesh->jump(in, tmp_npf0.dat);
  mesh->avg(tmp_grad0.dat, tmp_npf1.dat);
  mesh->avg(tmp_grad1.dat, tmp_npf2.dat);
  mesh->avg(tmp_grad2.dat, tmp_npf3.dat);
  timer->endTimer("PoissonMatrixFreeMult3D - mult faces");

  timer->startTimer("PoissonMatrixFreeMult3D - mult bfaces");
  if(mesh->bface2cells) {
    op_par_loop(pmf_3d_mult_avg_jump, "pmf_3d_mult_avg_jump", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mat_free_bcs, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_grad0.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_grad1.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_grad2.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf0.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf1.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf2.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf3.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonMatrixFreeMult3D - mult bfaces");

  timer->startTimer("PoissonMatrixFreeMult3D - finish flux");
  op_par_loop(pmf_3d_mult_flux, "pmf_3d_mult_flux", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->nx_c, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny_c, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz_c, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ_c, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_tau_c, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf0.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_npf1.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_npf2.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_npf3.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));
  timer->endTimer("PoissonMatrixFreeMult3D - finish flux");

  timer->startTimer("PoissonMatrixFreeMult3D - cells");
  timer->startTimer("PoissonMatrixFreeMult3D - mult cells MM");
  mesh->mass(tmp_grad0.dat);
  mesh->mass(tmp_grad1.dat);
  mesh->mass(tmp_grad2.dat);
  timer->endTimer("PoissonMatrixFreeMult3D - mult cells MM");
  timer->startTimer("PoissonMatrixFreeMult3D - mult cells Emat");
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf1.dat, 1.0, tmp_grad0.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf2.dat, 1.0, tmp_grad1.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf3.dat, 1.0, tmp_grad2.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf0.dat, 0.0, out);
  timer->endTimer("PoissonMatrixFreeMult3D - mult cells Emat");
  timer->startTimer("PoissonMatrixFreeMult3D - mult cells geof");
  op_par_loop(pmf_3d_mult_cells_geof, "pmf_3d_mult_cells_geof", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 10, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_grad0.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_grad1.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_grad2.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("PoissonMatrixFreeMult3D - mult cells geof");
  timer->startTimer("PoissonMatrixFreeMult3D - mult cells D");
  op2_gemv(mesh, true, 1.0, DGConstants::DR, tmp_grad0.dat, 1.0, out);
  op2_gemv(mesh, true, 1.0, DGConstants::DS, tmp_grad1.dat, 1.0, out);
  op2_gemv(mesh, true, 1.0, DGConstants::DT, tmp_grad2.dat, 1.0, out);
  timer->endTimer("PoissonMatrixFreeMult3D - mult cells D");
  timer->endTimer("PoissonMatrixFreeMult3D - cells");

  dg_dat_pool->releaseTempDatCells(tmp_npf0);
  dg_dat_pool->releaseTempDatCells(tmp_npf1);
  dg_dat_pool->releaseTempDatCells(tmp_npf2);
  dg_dat_pool->releaseTempDatCells(tmp_npf3);
  dg_dat_pool->releaseTempDatCells(tmp_grad0);
  dg_dat_pool->releaseTempDatCells(tmp_grad1);
  dg_dat_pool->releaseTempDatCells(tmp_grad2);
  timer->endTimer("PoissonMatrixFreeMult3D - mult");
}
