#include "matrices/3d/poisson_matrix_free_3d.h"

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
  op_arg arg16,
  op_arg arg17,
  op_arg arg18,
  op_arg arg19);

PoissonMatrixFree3D::PoissonMatrixFree3D(DGMesh3D *m) {
  mesh = m;
  _mesh = m;

  DG_FP *tmp_np  = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *tmp_npf_data = (DG_FP *)calloc(4 * DG_NPF * mesh->cells->size, sizeof(DG_FP));

  in_grad[0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_in_0");
  in_grad[1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_in_1");
  in_grad[2] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_in_2");
  l[0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_l0");
  l[1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_l1");
  l[2] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_l2");
  tmp_npf[0] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, tmp_npf_data, "poisson_matrix_free_tmp_npf0");
  tmp_npf[1] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, tmp_npf_data, "poisson_matrix_free_tmp_npf1");
  tmp_npf[2] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, tmp_npf_data, "poisson_matrix_free_tmp_npf2");
  tmp_npf[3] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, tmp_npf_data, "poisson_matrix_free_tmp_npf3");

  free(tmp_npf_data);
  free(tmp_np);
}

void PoissonMatrixFree3D::apply_bc(op_dat rhs, op_dat bc) {
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
}

/*
void PoissonMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonMatrixFree3D - mult");
  timer->startTimer("PoissonMatrixFree3D - mult grad");
  mesh->grad(in, in_grad[0], in_grad[1], in_grad[2]);
  timer->endTimer("PoissonMatrixFree3D - mult grad");

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

  timer->endTimer("PoissonMatrixFree3D - mult");
}
*/

void PoissonMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonMatrixFree3D - mult");
  timer->startTimer("PoissonMatrixFree3D - mult grad");
  mesh->grad(in, in_grad[0], in_grad[1], in_grad[2]);
  timer->endTimer("PoissonMatrixFree3D - mult grad");

  op_par_loop(zero_npf_1, "zero_npf_1", _mesh->cells,
              op_arg_dat(tmp_npf[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(zero_npf_3, "zero_npf_3", _mesh->cells,
              op_arg_dat(tmp_npf[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_npf[3], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("PoissonMatrixFree3D - mult faces flux");
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
  timer->endTimer("PoissonMatrixFree3D - mult faces flux");

  timer->startTimer("PoissonMatrixFree3D - mult faces bflux");
  if(mesh->bflux2cells) {
    op_par_loop(pmf_3d_mult_faces_bflux, "pmf_3d_mult_faces_bflux", mesh->fluxes,
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
                op_arg_dat(tmp_npf[0], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
                op_arg_dat(tmp_npf[1], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
                op_arg_dat(tmp_npf[2], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
                op_arg_dat(tmp_npf[3], 0, mesh->bflux2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("PoissonMatrixFree3D - mult faces bflux");
/*
  timer->startTimer("PoissonMatrixFree3D - mult faces");
  op_par_loop(pmf_3d_mult_faces, "pmf_3d_mult_faces", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(in, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[0], -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(tmp_npf[1], -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(tmp_npf[2], -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(tmp_npf[3], -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  timer->endTimer("PoissonMatrixFree3D - mult faces");
*/

  timer->startTimer("PoissonMatrixFree3D - mult bfaces");
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
  timer->endTimer("PoissonMatrixFree3D - mult bfaces");

  timer->startTimer("PoissonMatrixFree3D - mult Emat");
  #ifdef OP2_DG_CUDA
  custom_kernel_pmf_3d_mult_cells_emat("pmf_3d_mult_cells_emat", _mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::EMAT), DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[3], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(l[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(l[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(l[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(out,  -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  #else
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[0], 0.0, l[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[1], 0.0, l[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[2], 0.0, l[2]);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[3], 0.0, out);
/*
  op_par_loop(pmf_3d_mult_cells_emat, "pmf_3d_mult_cells_emat", _mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::EMAT), DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_npf[3], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(l[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(l[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(l[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(out,  -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
*/
  #endif
  timer->endTimer("PoissonMatrixFree3D - mult Emat");

  timer->startTimer("PoissonMatrixFree3D - mult MM");
  mesh->mass(in_grad[0]);
  mesh->mass(in_grad[1]);
  mesh->mass(in_grad[2]);
  timer->endTimer("PoissonMatrixFree3D - mult MM");

  timer->startTimer("PoissonMatrixFree3D - mult cells");
  #ifdef OP2_DG_CUDA
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
              op_arg_dat(l[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  #else
  op_par_loop(pmf_3d_mult_cells, "pmf_3d_mult_cells", _mesh->cells,
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
              op_arg_dat(l[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  #endif
  timer->endTimer("PoissonMatrixFree3D - mult cells");
  timer->endTimer("PoissonMatrixFree3D - mult");
}
