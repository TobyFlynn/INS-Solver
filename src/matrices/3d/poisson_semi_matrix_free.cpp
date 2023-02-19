#include "matrices/3d/poisson_semi_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

PoissonSemiMatrixFree3D::PoissonSemiMatrixFree3D(DGMesh3D *m) {
  mesh = m;
  _mesh = m;

  DG_FP *data_t0 = (DG_FP *)calloc(DG_NP * DG_NP * mesh->cells->size, sizeof(DG_FP));
  op1 = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, data_t0, "poisson_matrix_op1");
  free(data_t0);

  DG_FP *data_t1 = (DG_FP *)calloc(DG_NP * DG_NP * mesh->faces->size, sizeof(DG_FP));
  op2[0] = op_decl_dat(mesh->faces, DG_NP * DG_NP, DG_FP_STR, data_t1, "poisson_matrix_op20");
  op2[1] = op_decl_dat(mesh->faces, DG_NP * DG_NP, DG_FP_STR, data_t1, "poisson_matrix_op21");
  free(data_t1);

  DG_FP *data_t2 = (DG_FP *)calloc(DG_NP * DG_NPF * mesh->bfaces->size, sizeof(DG_FP));
  opbc = op_decl_dat(mesh->bfaces, DG_NP * DG_NPF, DG_FP_STR, data_t2, "poisson_matrix_opbc");
  free(data_t2);

  int *data_t3 = (int *)calloc(mesh->cells->size, sizeof(int));
  glb_ind = op_decl_dat(mesh->cells, 1, "int", data_t3, "poisson_matrix_glb_ind");
  free(data_t3);

  int *data_t4 = (int *)calloc(mesh->faces->size, sizeof(int));
  glb_indL = op_decl_dat(mesh->faces, 1, "int", data_t4, "poisson_matrix_glb_indL");
  glb_indR = op_decl_dat(mesh->faces, 1, "int", data_t4, "poisson_matrix_glb_indR");
  orderL  = op_decl_dat(mesh->faces, 1, "int", data_t4, "poisson_orderL");
  orderR  = op_decl_dat(mesh->faces, 1, "int", data_t4, "poisson_orderR");
  free(data_t4);

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

void PoissonSemiMatrixFree3D::calc_mat() {
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  petscMatResetRequired = true;
}

void PoissonSemiMatrixFree3D::apply_bc(op_dat rhs, op_dat bc) {
  if(mesh->bface2cells) {
    op_par_loop(poisson_matrix_3d_apply_bc, "poisson_matrix_3d_apply_bc", mesh->bfaces,
                op_arg_dat(opbc, -1, OP_ID, DG_NPF * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc,   -1, OP_ID, DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(rhs,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
}

void PoissonSemiMatrixFree3D::calc_op1() {
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
}

void PoissonSemiMatrixFree3D::calc_op2() {
  // TODO full p-adaptivity
  op_par_loop(poisson_matrix_3d_op2, "poisson_matrix_3d_op2", mesh->faces,
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
              op_arg_dat(op1, 1, mesh->face2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE));
}

void PoissonSemiMatrixFree3D::calc_opbc() {
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

    op_par_loop(poisson_matrix_3d_opbc, "poisson_matrix_3d_opbc", mesh->bfaces,
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
                op_arg_dat(opbc, -1, OP_ID, DG_NPF * DG_NP, DG_FP_STR, OP_WRITE));
  }
}

void PoissonSemiMatrixFree3D::calc_glb_ind() {
  set_glb_ind();
  op_par_loop(copy_to_edges, "copy_to_edges", mesh->faces,
              op_arg_dat(glb_ind, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));

  op_par_loop(copy_to_edges, "copy_to_edges", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_WRITE));
}

void PoissonSemiMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonSMF - Mult");
  timer->startTimer("PoissonSMF - Mult - Grad");
  mesh->grad(in, in_grad[0], in_grad[1], in_grad[2]);
  timer->endTimer("PoissonSMF - Mult - Grad");

  op_par_loop(zero_np_1, "zero_np_1", _mesh->cells,
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  op_par_loop(zero_np_3, "zero_np_3", _mesh->cells,
              op_arg_dat(l[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(l[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(l[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

op_par_loop(zero_npf_1, "zero_npf_1", _mesh->cells,
            op_arg_dat(tmp_npf[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

op_par_loop(zero_npf_3, "zero_npf_3", _mesh->cells,
            op_arg_dat(tmp_npf[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
            op_arg_dat(tmp_npf[2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
            op_arg_dat(tmp_npf[3], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("PoissonSMF - Mult - Faces");
  op_par_loop(pmf_3d_mult_faces_0, "pmf_3d_mult_faces_0", mesh->faces,
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
  timer->endTimer("PoissonSMF - Mult - Faces");
  // op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[0], 0.0, out);

  timer->startTimer("PoissonSMF - Mult - Emat");
  // op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[0], 0.0, l[0]);
  // op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[1], 0.0, l[1]);
  // op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[2], 0.0, l[2]);
  // op2_gemv(mesh, false, 1.0, DGConstants::EMAT, tmp_npf[3], 0.0, out);

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

  timer->endTimer("PoissonSMF - Mult - Emat");
  timer->startTimer("PoissonSMF - Mult - Cells");
  op_par_loop(pmf_3d_mult_cells, "pmf_3d_mult_cells", _mesh->cells,
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
              op_arg_dat(l[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("PoissonSMF - Mult - Cells");
  timer->endTimer("PoissonSMF - Mult");
/*
  timer->startTimer("PoissonMF - Mult");
  timer->startTimer("PoissonMF - Mult - Cells");
  op_par_loop(poisson_matrix_free_3d_mult_cells, "poisson_matrix_free_3d_mult_cells", _mesh->cells,
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
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonMF - Mult - Cells");
  timer->startTimer("PoissonMF - Mult - Faces");

  op_par_loop(poisson_matrix_free_3d_mult_faces, "poisson_matrix_free_3d_mult_faces", mesh->faces,
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
              op_arg_dat(in,  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_INC));

  timer->endTimer("PoissonMF - Mult - Faces");
  timer->endTimer("PoissonMF - Mult");
  */
}
