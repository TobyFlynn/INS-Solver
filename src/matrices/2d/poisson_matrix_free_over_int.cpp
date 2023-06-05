#include "matrices/2d/poisson_matrix_free_over_int_2d.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"

extern Timing *timer;
extern DGConstants *constants;

PoissonMatrixFreeOverInt2D::PoissonMatrixFreeOverInt2D(DGMesh2D *m) {
  mesh = m;
  _mesh = m;
  petscMatInit = false;

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

void PoissonMatrixFreeOverInt2D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonMatrixFreeOverInt2D - mult");
  timer->startTimer("PoissonMatrixFreeOverInt2D - mult grad");
  mesh->grad(in, in_grad[0], in_grad[1]);
  timer->endTimer("PoissonMatrixFreeOverInt2D - mult grad");

  timer->startTimer("PoissonMatrixFreeOverInt2D - mult faces");
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, in, 0.0, gIn);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, in_grad[0], 0.0, gIn_grad[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, in_grad[1], 0.0, gIn_grad[1]);

  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(g_tmp[0], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(g_tmp[1], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(g_tmp[2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

  op_par_loop(pmf_2d_mult_faces_over_int, "pmf_2d_mult_faces_over_int", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale_c, -2, mesh->face2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(gIn, -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gIn_grad[0], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(gIn_grad[1], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(g_tmp[0], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC),
              op_arg_dat(g_tmp[1], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC),
              op_arg_dat(g_tmp[2], -2, mesh->face2cells, DG_G_NP, DG_FP_STR, OP_INC));

  if(mesh->bface2cells) {
    op_par_loop(pmf_2d_mult_bfaces_over_int, "pmf_2d_mult_bfaces_over_int", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->fscale_c, 0, mesh->bface2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
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
  timer->endTimer("PoissonMatrixFreeOverInt2D - mult faces");

  timer->startTimer("PoissonMatrixFreeOverInt2D - mult MM");
  mesh->mass(in_grad[0]);
  mesh->mass(in_grad[1]);
  timer->endTimer("PoissonMatrixFreeOverInt2D - mult MM");

  timer->startTimer("PoissonMatrixFreeOverInt2D - mult cells");
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
  timer->endTimer("PoissonMatrixFreeOverInt2D - mult cells");
  timer->endTimer("PoissonMatrixFreeOverInt2D - mult");
}

void PoissonMatrixFreeOverInt2D::apply_bc(op_dat rhs, op_dat bc) {
  timer->startTimer("PoissonMatrixFreeOverInt2D - apply_bc");
  if(mesh->bface2cells) {
    op_par_loop(pmf_2d_apply_bc_over_int, "pmf_2d_apply_bc_over_int", mesh->bfaces,
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
                op_arg_dat(mesh->fscale_c, 0, mesh->bface2cells, 3 * DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(bc,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rhs, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonMatrixFreeOverInt2D - apply_bc");
}
