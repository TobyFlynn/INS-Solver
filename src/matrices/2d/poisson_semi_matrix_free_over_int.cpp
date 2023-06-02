#include "matrices/2d/poisson_semi_matrix_free_over_int_2d.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"

extern Timing *timer;
extern DGConstants *constants;

PoissonSemiMatrixFreeOverInt2D::PoissonSemiMatrixFreeOverInt2D(DGMesh2D *m) {
  mesh = m;
  _mesh = m;
  petscMatInit = false;

  DG_FP *tmp_np_np_c = (DG_FP *)calloc(DG_NP * DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *tmp_np_np_e = (DG_FP *)calloc(DG_NP * DG_NP * mesh->faces->size, sizeof(DG_FP));
  DG_FP *tmp_gf_np_be = (DG_FP *)calloc(DG_GF_NP * DG_NP * mesh->bfaces->size, sizeof(DG_FP));
  DG_FP *tmp_1 = (DG_FP *)calloc(mesh->cells->size, sizeof(DG_FP));
  int *tmp_1_int_c = (int *)calloc(mesh->cells->size, sizeof(int));
  int *tmp_1_int_e = (int *)calloc(mesh->faces->size, sizeof(int));
  int *tmp_1_int_be = (int *)calloc(mesh->bfaces->size, sizeof(int));

  op1      = op_decl_dat(mesh->cells, DG_NP * DG_NP, DG_FP_STR, tmp_np_np_c, "poisson_op1");
  h        = op_decl_dat(mesh->cells, 1, DG_FP_STR, tmp_1, "poisson_h");

  glb_ind   = op_decl_dat(mesh->cells, 1, "int", tmp_1_int_c, "poisson_glb_ind");
  glb_indL  = op_decl_dat(mesh->faces, 1, "int", tmp_1_int_e, "poisson_glb_indL");
  glb_indR  = op_decl_dat(mesh->faces, 1, "int", tmp_1_int_e, "poisson_glb_indR");
  glb_indBC = op_decl_dat(mesh->bfaces, 1, "int", tmp_1_int_be, "poisson_glb_indBC");

  orderL  = op_decl_dat(mesh->faces, 1, "int", tmp_1_int_e, "poisson_orderL");
  orderR  = op_decl_dat(mesh->faces, 1, "int", tmp_1_int_e, "poisson_orderR");
  orderBC = op_decl_dat(mesh->bfaces, 1, "int", tmp_1_int_be, "poisson_orderBC");

  free(tmp_1_int_be);
  free(tmp_1_int_e);
  free(tmp_1_int_c);
  free(tmp_1);
  free(tmp_gf_np_be);
  free(tmp_np_np_e);
  free(tmp_np_np_c);

  DG_FP *tmp_np  = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *tmp_g_np_data = (DG_FP *)calloc(DG_G_NP * mesh->cells->size, sizeof(DG_FP));

  in_grad[0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_in_0");
  in_grad[1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_in_1");
  gIn = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np_data, "poisson_matrix_free_gIn");
  gIn_grad[0] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np_data, "poisson_matrix_free_gIn_grad0");
  gIn_grad[1] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np_data, "poisson_matrix_free_gIn_grad1");
  g_tmp[0] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np_data, "poisson_matrix_free_g_tmp0");
  g_tmp[1] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np_data, "poisson_matrix_free_g_tmp1");
  g_tmp[2] = op_decl_dat(mesh->cells, DG_G_NP, DG_FP_STR, tmp_g_np_data, "poisson_matrix_free_g_tmp2");
  l[0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_l0");
  l[1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "poisson_matrix_free_l1");
  free(tmp_g_np_data);
  free(tmp_np);
}

void PoissonSemiMatrixFreeOverInt2D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - mult");
  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - mult grad");
  mesh->grad(in, in_grad[0], in_grad[1]);
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - mult grad");

  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - mult faces");
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, in, 0.0, gIn);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, in_grad[0], 0.0, gIn_grad[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, in_grad[1], 0.0, gIn_grad[1]);

  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(g_tmp[0], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(g_tmp[1], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(g_tmp[2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

  op_par_loop(pmf_2d_mult_faces, "pmf_2d_mult_faces", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -2, mesh->face2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(gIn, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gIn_grad[0], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gIn_grad[1], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(g_tmp[0], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC),
              op_arg_dat(g_tmp[1], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC),
              op_arg_dat(g_tmp[2], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC));

  if(mesh->bface2cells) {
    op_par_loop(pmf_2d_mult_bfaces, "pmf_2d_mult_bfaces", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->fscale, 0, mesh->bface2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(gIn, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gIn_grad[0], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gIn_grad[1], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(g_tmp[0], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC),
                op_arg_dat(g_tmp[1], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC),
                op_arg_dat(g_tmp[2], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, true, 1.0, DGConstants::GAUSS_INTERP, g_tmp[0], 0.0, out);
  op2_gemv(mesh, true, 1.0, DGConstants::GAUSS_INTERP, g_tmp[1], 0.0, l[0]);
  op2_gemv(mesh, true, 1.0, DGConstants::GAUSS_INTERP, g_tmp[2], 0.0, l[1]);
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - mult faces");

  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - mult MM");
  mesh->mass(in_grad[0]);
  mesh->mass(in_grad[1]);
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - mult MM");

  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - mult cells");
  op_par_loop(pmf_2d_mult_cells, "pmf_2d_mult_cells", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(l[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in_grad[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - mult cells");
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - mult");
}

void PoissonSemiMatrixFreeOverInt2D::calc_mat_partial() {
  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - calc_mat_partial");
  op_par_loop(poisson_h, "poisson_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(h, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

  // calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  petscMatResetRequired = true;
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - calc_mat_partial");
}

void PoissonSemiMatrixFreeOverInt2D::calc_glb_ind() {
  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - calc_glb_ind");
  set_glb_ind();
  op_par_loop(copy_to_edges, "copy_to_edges", mesh->faces,
              op_arg_dat(glb_ind, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));
  if(mesh->bface2cells) {
    op_par_loop(copy_to_bedges, "copy_to_bedges", mesh->bfaces,
                op_arg_dat(glb_ind, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(glb_indBC, -1, OP_ID, 1, "int", OP_WRITE));
  }

  op_par_loop(copy_to_edges, "copy_to_edges", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_WRITE));
  if(mesh->bface2cells) {
    op_par_loop(copy_to_bedges, "copy_to_bedges", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(orderBC, -1, OP_ID, 1, "int", OP_WRITE));
  }
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - calc_glb_ind");
}

void PoissonSemiMatrixFreeOverInt2D::apply_bc(op_dat rhs, op_dat bc) {
  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - apply_bc");
  if(mesh->bface2cells) {
    op_par_loop(pmf_2d_apply_bc, "pmf_2d_apply_bc", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP0), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP1), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP2), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->fscale, 0, mesh->bface2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(bc,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rhs, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - apply_bc");
}

void PoissonSemiMatrixFreeOverInt2D::calc_op1() {
  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - calc_op1");
  op_par_loop(poisson_matrix_2d_op1, "poisson_matrix_2d_op1", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - calc_op1");
}

void PoissonSemiMatrixFreeOverInt2D::calc_op2() {
  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - calc_op2");
  op_par_loop(poisson_matrix_2d_op2_partial, "poisson_matrix_2d_op2_partial", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP0), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP1), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP2), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->x, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -2, mesh->face2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(op1, 0, mesh->face2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC),
              op_arg_dat(op1, 1, mesh->face2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC));
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - calc_op2");
}

void PoissonSemiMatrixFreeOverInt2D::calc_opbc() {
  timer->startTimer("PoissonSemiMatrixFreeOverInt2D - calc_opbc");
  if(mesh->bface2cells) {
    op_par_loop(poisson_matrix_2d_bop, "poisson_matrix_2d_bop", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DR), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DS), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP0), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP1), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP2), DG_ORDER * DG_GF_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->fscale, 0, mesh->bface2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(op1, 0, mesh->bface2cells, DG_NP * DG_NP, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonSemiMatrixFreeOverInt2D - calc_opbc");
}
