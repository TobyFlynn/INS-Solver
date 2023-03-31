#include "matrices/3d/factor_poisson_matrix_free_mult_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

void custom_kernel_pmf_3d_mult_cells_emat(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9);

void custom_kernel_pmf_3d_mult_cells(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16);

void custom_kernel_fpmf_grad_3d(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg argFactor,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16);

void custom_kernel_fpmf_3d_mult_faces_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg17,
  op_arg arg18,
  op_arg arg22,
  op_arg arg26,
  op_arg arg30,
  op_arg arg31,
  op_arg arg32,
  op_arg arg33);

FactorPoissonMatrixFreeMult3D::FactorPoissonMatrixFreeMult3D(DGMesh3D *m, bool alloc_tmp_dats) : PoissonMatrixFreeMult3D(m, alloc_tmp_dats) {
  DG_FP *tmp_4 = (DG_FP *)calloc(4 * mesh->fluxes->size, sizeof(DG_FP));
  mat_free_gtau = op_decl_dat(mesh->fluxes, 4, DG_FP_STR, tmp_4, "poisson_matrix_free_tau");
  free(tmp_4);
}

void FactorPoissonMatrixFreeMult3D::mat_free_set_factor(op_dat f) {
  mat_free_factor = f;

  timer->startTimer("FactorPoissonMatrixFreeMult3D - calc tau");
  op_par_loop(fpmf_3d_calc_tau, "fpmf_3d_calc_tau", mesh->fluxes,
              op_arg_dat(mesh->order, 0, mesh->flux2main_cell, 1, "int", OP_READ),
              op_arg_dat(mesh->fluxFaceNums, -1, OP_ID, 8, "int", OP_READ),
              op_arg_dat(mesh->fluxFmask,    -1, OP_ID, 4 * DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fluxFscale, -1, OP_ID, 8, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor, 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor, -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_gtau, -1, OP_ID, 4, DG_FP_STR, OP_WRITE));
  timer->endTimer("FactorPoissonMatrixFreeMult3D - calc tau");
}

void FactorPoissonMatrixFreeMult3D::mat_free_apply_bc(op_dat rhs, op_dat bc) {
  if(mesh->bface2cells) {
    op_par_loop(fpmf_3d_apply_bc, "fpmf_3d_apply_bc", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F0), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F1), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F2), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F3), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mat_free_bcs, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->rx, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sx, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->tx, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->ry, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sy, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->ty, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->rz, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sz, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->tz, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_factor,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(rhs, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
}

void FactorPoissonMatrixFreeMult3D::mat_free_mult(op_dat in, op_dat out) {
  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult");
  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult grad");
  #ifdef OP2_DG_CUDA
  custom_kernel_fpmf_grad_3d("fpmf_grad_3d", mesh->cells,
                       op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
                       op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                       op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                       op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                       op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                       op_arg_dat(mat_free_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->rx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->sx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->tx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->ry, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->sy, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->ty, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->rz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->sz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(mesh->tz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                       op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                       op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                       op_arg_dat(in_grad[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  #else
    op_par_loop(fpmf_grad_3d, "fpmf_grad_3d", mesh->cells,
                op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->rx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->tx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->ry, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sy, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->ty, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->rz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->tz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(in_grad[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  #endif
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult grad");

  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(tmp_npf[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(tmp_npf[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[3], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult faces flux");
  #ifdef OP2_DG_CUDA
  custom_kernel_fpmf_3d_mult_faces_flux("fpmf_3d_mult_faces_flux", mesh->fluxes,
              op_arg_dat(mesh->order, 0, mesh->flux2main_cell, 1, "int", OP_READ),
              op_arg_dat(mesh->fluxFaceNums, -1, OP_ID, 8, "int", OP_READ),
              op_arg_dat(mesh->fluxFmask,    -1, OP_ID, 4 * DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fluxNx, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxNy, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxNz, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxFscale, -1, OP_ID, 8, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxSJ, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_gtau, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(in, 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in, -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor, 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[0], 0, mesh->flux2main_cell, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[1], 0, mesh->flux2main_cell, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[2], 0, mesh->flux2main_cell, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[3], 0, mesh->flux2main_cell, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  #else
  op_par_loop(fpmf_3d_mult_faces_flux, "fpmf_3d_mult_faces_flux", mesh->fluxes,
              op_arg_dat(mesh->order, 0, mesh->flux2main_cell, 1, "int", OP_READ),
              op_arg_dat(mesh->fluxFaceNums, -1, OP_ID, 8, "int", OP_READ),
              op_arg_dat(mesh->fluxFmask,    -1, OP_ID, 4 * DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fluxNx, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxNy, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxNz, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxFscale, -1, OP_ID, 8, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxSJ, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_gtau, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(in, 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in, -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor, 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], 0, mesh->flux2main_cell, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -4, mesh->flux2neighbour_cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[0], 0, mesh->flux2main_cell, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[1], 0, mesh->flux2main_cell, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[2], 0, mesh->flux2main_cell, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[3], 0, mesh->flux2main_cell, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  #endif
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult faces flux");

  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult faces bflux");
  if(mesh->bflux2cells) {
    op_par_loop(fpmf_3d_mult_faces_bflux, "fpmf_3d_mult_faces_bflux", mesh->bfluxes,
                op_arg_dat(mesh->order,   0, mesh->bflux2cells, 1, "int", OP_READ),
                op_arg_dat(mesh->bfluxL, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->faceNum, 0, mesh->bflux2faces, 2, "int", OP_READ),
                op_arg_dat(mesh->fmaskL,  0, mesh->bflux2faces, DG_NPF, "int", OP_READ),
                op_arg_dat(mesh->fmaskR,  0, mesh->bflux2faces, DG_NPF, "int", OP_READ),
                op_arg_dat(mesh->nx, 0, mesh->bflux2faces, 2, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->ny, 0, mesh->bflux2faces, 2, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->nz, 0, mesh->bflux2faces, 2, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->fscale, 0, mesh->bflux2faces, 2, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sJ, 0, mesh->bflux2faces, 2, DG_FP_STR, OP_READ),
                op_arg_dat(in, -2, mesh->bflux2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_factor, -2, mesh->bflux2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[0], -2, mesh->bflux2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[1], -2, mesh->bflux2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[2], -2, mesh->bflux2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf[0], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[1], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[2], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[3], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult faces bflux");

  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult bfaces");
  if(mesh->bface2cells) {
    op_par_loop(fpmf_3d_mult_bfaces, "fpmf_3d_mult_bfaces", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(mat_free_bcs, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_factor, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf[0], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[1], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[2], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[3], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult bfaces");

  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult cells");
  #ifdef OP2_DG_CUDA
  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult cells MM");
  mesh->mass(in_grad[0]);
  mesh->mass(in_grad[1]);
  mesh->mass(in_grad[2]);
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult cells MM");
  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult cells Emat");
  custom_kernel_pmf_3d_mult_cells_emat("pmf_3d_mult_cells_emat", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::EMAT), DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[3], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(in_grad[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(out,  -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult cells Emat");
  timer->startTimer("FactorPoissonMatrixFreeMult3D - mult cells cells");
  custom_kernel_pmf_3d_mult_cells("pmf_3d_mult_cells", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ty, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult cells cells");
  #else
  op_par_loop(pmf_3d_mult_cells_merged, "pmf_3d_mult_cells_merged", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::EMAT), DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ty, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[3], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  #endif
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult cells");
  timer->endTimer("FactorPoissonMatrixFreeMult3D - mult");
}
