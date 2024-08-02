#include "slip_matrix/2d/viscous_matrix.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"
#include "dg_dat_pool.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;
extern DGDatPool *dg_dat_pool;

ViscousMatrix2D::ViscousMatrix2D(DGMesh2D *m) {
  mesh = m;

  mat_free_tau_c = op_decl_dat(mesh->cells, 3, DG_FP_STR, (DG_FP *)NULL, "mat_free_tau_c");
}

void ViscousMatrix2D::set_bc_types(op_dat u_bc_ty, op_dat v_bc_ty) {
  u_bc_types = u_bc_ty;
  v_bc_types = v_bc_ty;
}

void ViscousMatrix2D::apply_bc(op_dat u_rhs, op_dat v_rhs, op_dat u_bc, op_dat v_bc) {
  if(mesh->bface2cells) {
    op_par_loop(vmf_2d_apply_bc, "vmf_2d_apply_bc", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(u_bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(v_bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, 0, mesh->bface2cells, 5, DG_FP_STR, OP_READ),
                op_arg_dat(u_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(v_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(u_rhs, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC),
                op_arg_dat(v_rhs, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
}

void ViscousMatrix2D::set_factor(DG_FP f) {
  factor = f;
}

DG_FP ViscousMatrix2D::get_factor() {
  return factor;
}

void ViscousMatrix2D::mat_free_pre_compute_tau() {
  timer->startTimer("ViscousMatrix2D - calc tau");
  op_par_loop(vmf_2d_calc_tau_faces, "vmf_2d_calc_tau_faces", mesh->faces,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_tau_c, -2, mesh->face2cells, 3, DG_FP_STR, OP_WRITE));
  if(mesh->bface2cells) {
    op_par_loop(vmf_2d_calc_tau_bfaces, "vmf_2d_calc_tau_bfaces", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_tau_c, 0, mesh->bface2cells, 3, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("ViscousMatrix2D - calc tau");
}

void ViscousMatrix2D::mult(op_dat u_in, op_dat v_in, op_dat u_out, op_dat v_out) {
  mat_free_pre_compute_tau();
  timer->startTimer("ViscousMatrix2D - mult");
  DGTempDat u_x = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat u_y = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat v_x = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat v_y = dg_dat_pool->requestTempDatCells(DG_NP);

  mesh->grad(u_in, u_x.dat, u_y.dat);
  mesh->grad(v_in, v_x.dat, v_y.dat);

  DGTempDat u_jump  = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat u_x_avg = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat u_y_avg = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat v_jump  = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat v_x_avg = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat v_y_avg = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(u_jump.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(u_x_avg.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(u_y_avg.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  
  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(v_jump.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(v_x_avg.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(v_y_avg.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  
  mesh->jump(u_in, u_jump.dat);
  mesh->avg(u_x.dat, u_x_avg.dat);
  mesh->avg(u_y.dat, u_y_avg.dat);
  mesh->jump(v_in, v_jump.dat);
  mesh->avg(v_x.dat, v_x_avg.dat);
  mesh->avg(v_y.dat, v_y_avg.dat);

  if(mesh->bface2cells) {
    op_par_loop(vmf_2d_mult_avg_jump, "vmf_2d_mult_avg_jump", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(u_bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(v_bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(u_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u_x.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u_y.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v_x.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v_y.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u_jump.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(u_x_avg.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(u_y_avg.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(v_jump.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(v_x_avg.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(v_y_avg.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  op_par_loop(vmf_2d_mult_flux, "vmf_2d_mult_flux", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->nx_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_tau_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(u_jump.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(u_x_avg.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(u_y_avg.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(v_jump.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(v_x_avg.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(v_y_avg.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));

  mesh->mass(u_x.dat);
  mesh->mass(u_y.dat);
  mesh->mass(v_x.dat);
  mesh->mass(v_y.dat);

  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, u_x_avg.dat, 1.0, u_x.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, u_y_avg.dat, 1.0, u_y.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, u_jump.dat, 0.0, u_out);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, v_x_avg.dat, 1.0, v_x.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, v_y_avg.dat, 1.0, v_y.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, v_jump.dat, 0.0, v_out);

  op_par_loop(vmf_2d_mult_cells_geof, "vmf_2d_mult_cells_geof", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(u_x.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(u_y.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(v_x.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(v_y.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, true, 1.0, DGConstants::DR, u_x.dat, 1.0, u_out);
  op2_gemv(mesh, true, 1.0, DGConstants::DS, u_y.dat, 1.0, u_out);
  op2_gemv(mesh, true, 1.0, DGConstants::DR, v_x.dat, 1.0, v_out);
  op2_gemv(mesh, true, 1.0, DGConstants::DS, v_y.dat, 1.0, v_out);

  op_par_loop(vmf_2d_mult_mm_geof, "vmf_2d_mult_mm_geof", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_gbl(&factor,  1, DG_FP_STR, OP_READ),
              op_arg_dat(u_in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v_in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u_out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(v_out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  dg_dat_pool->releaseTempDatCells(u_x);
  dg_dat_pool->releaseTempDatCells(u_y);
  dg_dat_pool->releaseTempDatCells(v_x);
  dg_dat_pool->releaseTempDatCells(v_y);
  dg_dat_pool->releaseTempDatCells(u_jump);
  dg_dat_pool->releaseTempDatCells(u_x_avg);
  dg_dat_pool->releaseTempDatCells(u_y_avg);
  dg_dat_pool->releaseTempDatCells(v_jump);
  dg_dat_pool->releaseTempDatCells(v_x_avg);
  dg_dat_pool->releaseTempDatCells(v_y_avg);
  timer->endTimer("ViscousMatrix2D - mult");
}