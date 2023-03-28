#include "matrices/3d/poisson_semi_matrix_free_3d.h"

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

PoissonSemiMatrixFree3D::PoissonSemiMatrixFree3D(DGMesh3D *m, bool alloc_tmp_dats) {
  mesh = m;
  _mesh = m;

  DG_FP *data_t0 = (DG_FP *)calloc(DG_NP * DG_NP * mesh->cells->size, sizeof(DG_FP));
  op1 = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, data_t0, "poisson_matrix_op1");
  free(data_t0);

  if(alloc_tmp_dats) {
    DG_FP *tmp_np  = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
    DG_FP *tmp_npf_data = (DG_FP *)calloc(4 * DG_NPF * mesh->cells->size, sizeof(DG_FP));

    in_grad[0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_in_0");
    in_grad[1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_in_1");
    in_grad[2] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_in_2");
    tmp_npf[0] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, tmp_npf_data, "poisson_matrix_free_tmp_npf0");
    tmp_npf[1] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, tmp_npf_data, "poisson_matrix_free_tmp_npf1");
    tmp_npf[2] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, tmp_npf_data, "poisson_matrix_free_tmp_npf2");
    tmp_npf[3] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, tmp_npf_data, "poisson_matrix_free_tmp_npf3");

    free(tmp_npf_data);
    free(tmp_np);
  }

  DG_FP *tmp_4 = (DG_FP *)calloc(4 * mesh->fluxes->size, sizeof(DG_FP));
  gtau = op_decl_dat(mesh->fluxes, 4, DG_FP_STR, tmp_4, "poisson_matrix_free_tau");
  free(tmp_4);
}

void PoissonSemiMatrixFree3D::set_tmp_dats(op_dat np0, op_dat np1, op_dat np2,
                                           op_dat npf0, op_dat npf1,
                                           op_dat npf2, op_dat npf3) {
  in_grad[0] = np0;
  in_grad[1] = np1;
  in_grad[2] = np2;
  tmp_npf[0] = npf0;
  tmp_npf[1] = npf1;
  tmp_npf[2] = npf2;
  tmp_npf[3] = npf3;
}

void PoissonSemiMatrixFree3D::calc_mat_partial() {
  timer->startTimer("PoissonSemiMatrixFree3D - calc_mat_partial");
  calc_op1();
  calc_op2();
  calc_opbc();
  petscMatResetRequired = true;
  timer->endTimer("PoissonSemiMatrixFree3D - calc_mat_partial");
}

void PoissonSemiMatrixFree3D::apply_bc(op_dat rhs, op_dat bc) {
  timer->startTimer("PoissonSemiMatrixFree3D - apply_bc");
  if(mesh->bface2cells) {
    op_par_loop(pmf_3d_apply_bc, "pmf_3d_apply_bc", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F0), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F1), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F2), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F3), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
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
                op_arg_dat(bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(rhs, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonSemiMatrixFree3D - apply_bc");
}

/*
void PoissonSemiMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonSemiMatrixFree3D - mult");
  timer->startTimer("PoissonSemiMatrixFree3D - mult grad");
  mesh->grad(in, in_grad[0], in_grad[1], in_grad[2]);
  timer->endTimer("PoissonSemiMatrixFree3D - mult grad");

  op_par_loop(pmf_3d_mult_complete_flux, "pmf_3d_mult_complete_flux", mesh->fluxes,
              op_arg_dat(mesh->order, 0, mesh->flux2cells, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::EMAT), DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxFaceNums, -1, OP_ID, 8, "int", OP_READ),
              op_arg_dat(mesh->fluxFmask,    -1, OP_ID, 4 * DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fluxNx, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxNy, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxNz, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxFscale, -1, OP_ID, 8, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxSJ, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rx, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sx, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tx, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ry, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sy, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ty, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rz, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sz, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tz, 0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J,  0, mesh->flux2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(in, -5, mesh->flux2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -5, mesh->flux2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -5, mesh->flux2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -5, mesh->flux2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, 0, mesh->flux2cells, DG_NP, DG_FP_STR, OP_WRITE));

  // TODO boundary conditions

  timer->endTimer("PoissonSemiMatrixFree3D - mult");
}
*/

void PoissonSemiMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonSemiMatrixFree3D - mult");
  timer->startTimer("PoissonSemiMatrixFree3D - mult grad");
  mesh->grad(in, in_grad[0], in_grad[1], in_grad[2]);
  timer->endTimer("PoissonSemiMatrixFree3D - mult grad");

  op_par_loop(zero_npf_1, "zero_npf_1", _mesh->cells,
              op_arg_dat(tmp_npf[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(zero_npf_3, "zero_npf_3", _mesh->cells,
              op_arg_dat(tmp_npf[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[3], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("PoissonSemiMatrixFree3D - mult faces flux");
  op_par_loop(pmf_3d_mult_faces_flux, "pmf_3d_mult_faces_flux", mesh->fluxes,
              op_arg_dat(mesh->order, 0, mesh->flux2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->fluxFaceNums, -1, OP_ID, 8, "int", OP_READ),
              op_arg_dat(mesh->fluxFmask,    -1, OP_ID, 4 * DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fluxNx, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxNy, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxNz, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxFscale, -1, OP_ID, 8, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fluxSJ, -1, OP_ID, 4, DG_FP_STR, OP_READ),
              op_arg_dat(in, -5, mesh->flux2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -5, mesh->flux2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -5, mesh->flux2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -5, mesh->flux2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[0], 0, mesh->flux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[1], 0, mesh->flux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[2], 0, mesh->flux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[3], 0, mesh->flux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonSemiMatrixFree3D - mult faces flux");

  timer->startTimer("PoissonSemiMatrixFree3D - mult faces bflux");
  if(mesh->bflux2cells) {
    op_par_loop(pmf_3d_mult_faces_bflux, "pmf_3d_mult_faces_bflux", mesh->bfluxes,
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
                op_arg_dat(in_grad[0], -2, mesh->bflux2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[1], -2, mesh->bflux2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[2], -2, mesh->bflux2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf[0], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[1], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[2], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[3], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonSemiMatrixFree3D - mult faces bflux");

  timer->startTimer("PoissonSemiMatrixFree3D - mult bfaces");
  if(mesh->bface2cells) {
    op_par_loop(pmf_3d_mult_bfaces, "pmf_3d_mult_bfaces", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(in_grad[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf[0], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[1], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[2], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_npf[3], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonSemiMatrixFree3D - mult bfaces");

  timer->startTimer("PoissonSemiMatrixFree3D - mult cells");
  #ifdef OP2_DG_CUDA
  timer->startTimer("PoissonMatrixFree3D - mult cells MM");
  mesh->mass(in_grad[0]);
  mesh->mass(in_grad[1]);
  mesh->mass(in_grad[2]);
  timer->endTimer("PoissonMatrixFree3D - mult cells MM");
  timer->startTimer("PoissonMatrixFree3D - mult cells Emat");
  custom_kernel_pmf_3d_mult_cells_emat("pmf_3d_mult_cells_emat", _mesh->cells,
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
  timer->endTimer("PoissonMatrixFree3D - mult cells Emat");
  timer->startTimer("PoissonMatrixFree3D - mult cells cells");
  custom_kernel_pmf_3d_mult_cells("pmf_3d_mult_cells", _mesh->cells,
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
  timer->endTimer("PoissonMatrixFree3D - mult cells cells");
  #else
  op_par_loop(pmf_3d_mult_cells_merged, "pmf_3d_mult_cells_merged", _mesh->cells,
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
  timer->endTimer("PoissonSemiMatrixFree3D - mult cells");
  timer->endTimer("PoissonSemiMatrixFree3D - mult");
}

void PoissonSemiMatrixFree3D::calc_op1() {
  timer->startTimer("PoissonSemiMatrixFree3D - calc_op1");
  op_par_loop(poisson_matrix_3d_op1, "poisson_matrix_3d_op1", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ty, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J,  -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonSemiMatrixFree3D - calc_op1");
}

void PoissonSemiMatrixFree3D::calc_op2() {
  timer->startTimer("PoissonSemiMatrixFree3D - calc_op2");
  // TODO full p-adaptivity
  op_par_loop(poisson_matrix_3d_op2_partial, "poisson_matrix_3d_op2_partial", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F0), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F1), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F2), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F3), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rx, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sx, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tx, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ry, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sy, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ty, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rz, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sz, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->tz, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(op1, 0, mesh->face2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC),
              op_arg_dat(op1, 1, mesh->face2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC));
  timer->endTimer("PoissonSemiMatrixFree3D - calc_op2");
}

void PoissonSemiMatrixFree3D::calc_opbc() {
  timer->startTimer("PoissonSemiMatrixFree3D - calc_opbc");
  if(mesh->bface2cells) {
    op_par_loop(poisson_matrix_3d_bop, "poisson_matrix_3d_bop", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F0), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F1), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F2), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F3), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
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
                op_arg_dat(mesh->tz, 0, mesh->bface2cells, 1, DG_FP_STR, OP_READ)
                op_arg_dat(op1, 0, mesh->bface2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonSemiMatrixFree3D - calc_opbc");
}
