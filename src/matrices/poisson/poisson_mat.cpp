#include "matrices/2d/poisson_matrix.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"

extern Timing *timer;
extern DGConstants *constants;

PoissonMatrix2D::PoissonMatrix2D(DGMesh2D *m) {
  mesh = m;

  double *tmp_np_np_c = (double *)calloc(DG_NP * DG_NP * mesh->cells->size, sizeof(double));
  double *tmp_np_np_e = (double *)calloc(DG_NP * DG_NP * mesh->faces->size, sizeof(double));
  double *tmp_gf_np_be = (double *)calloc(DG_GF_NP * DG_NP * mesh->bfaces->size, sizeof(double));
  double *tmp_1 = (double *)calloc(mesh->cells->size, sizeof(double));
  double *tmp_g_np = (double *)calloc(DG_G_NP * mesh->cells->size, sizeof(double));
  double *tmp_cub_np = (double *)calloc(DG_CUB_NP * mesh->cells->size, sizeof(double));
  int *tmp_1_int_c = (int *)calloc(mesh->cells->size, sizeof(int));
  int *tmp_1_int_e = (int *)calloc(mesh->faces->size, sizeof(int));
  int *tmp_1_int_be = (int *)calloc(mesh->bfaces->size, sizeof(int));

  op1      = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", tmp_np_np_c, "poisson_op1");
  op2[0]   = op_decl_dat(mesh->faces, DG_NP * DG_NP, "double", tmp_np_np_e, "poisson_op20");
  op2[1]   = op_decl_dat(mesh->faces, DG_NP * DG_NP, "double", tmp_np_np_e, "poisson_op21");
  op_bc    = op_decl_dat(mesh->bfaces, DG_GF_NP * DG_NP, "double", tmp_gf_np_be, "poisson_op_bc");
  h        = op_decl_dat(mesh->cells, 1, "double", tmp_1, "poisson_h");
  gFactor  = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_g_np, "poisson_gFactor");
  cFactor  = op_decl_dat(mesh->cells, DG_CUB_NP, "double", tmp_cub_np, "poisson_cFactor");

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
  free(tmp_cub_np);
  free(tmp_g_np);
  free(tmp_1);
  free(tmp_gf_np_be);
  free(tmp_np_np_e);
  free(tmp_np_np_c);
}

void PoissonMatrix2D::init() {
  op_par_loop(poisson_h, "poisson_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(h, -1, OP_ID, 1, "double", OP_WRITE));
}

void PoissonMatrix2D::calc_mat(op_dat fact) {
  timer->startTimer("PoissonMat - calc mat");
  factor = fact;
  calc_cub_sub_mat();
  calc_gauss_sub_mat();
  timer->endTimer("PoissonMat - calc mat");
}

void PoissonMatrix2D::calc_mat_mm(op_dat fact, op_dat mmFact) {
  timer->startTimer("PoissonMat - calc mat mm");
  factor = fact;
  mmFactor = mmFact;
  calc_cub_sub_mat();
  calc_gauss_sub_mat();
  calc_mm_mat();
  timer->endTimer("PoissonMat - calc mat mm");
}

void PoissonMatrix2D::update_glb_ind() {
  unknowns = get_local_unknowns();
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
}

void PoissonMatrix2D::mult(op_dat in, op_dat out) {
  op_par_loop(poisson_cells, "poisson_cells", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_WRITE));

  op_par_loop(poisson_edges, "poisson_edges", mesh->faces,
              op_arg_dat(mesh->order, 0, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(in,      0, mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out,     0, mesh->face2cells, DG_NP, "double", OP_INC),
              op_arg_dat(mesh->order, 1, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(in,      1, mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out,     1, mesh->face2cells, DG_NP, "double", OP_INC));
}

void PoissonMatrix2D::multJacobi(op_dat in, op_dat out) {
  op_par_loop(poisson_cells, "poisson_cells", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_WRITE));

  op_par_loop(poisson_edges, "poisson_edges", mesh->faces,
              op_arg_dat(mesh->order, 0, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(in,      0, mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out,     0, mesh->face2cells, DG_NP, "double", OP_INC),
              op_arg_dat(mesh->order, 1, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(in,      1, mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out,     1, mesh->face2cells, DG_NP, "double", OP_INC));

  op_par_loop(poisson_jacobi, "poisson_jacobi", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_RW));
}

void PoissonMatrix2D::transpose() {
  op_par_loop(transpose_cells, "transpose_cells", mesh->cells,
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_RW));

  op_par_loop(transpose_edges, "transpose_edges", mesh->faces,
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_RW),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_RW));
}

void PoissonMatrix2D::setDirichletBCs(int *d) {
  dirichlet[0] = d[0];
  dirichlet[1] = d[1];
  dirichlet[2] = d[2];
}

void PoissonMatrix2D::setNeumannBCs(int *n) {
  neumann[0] = n[0];
  neumann[1] = n[1];
  neumann[2] = n[2];
}

void PoissonMatrix2D::calc_cub_sub_mat() {
  // Initialise geometric factors for calcuating grad matrix
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDR, mesh->x, 0.0, mesh->cubature->op_tmp[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDS, mesh->x, 0.0, mesh->cubature->op_tmp[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDR, mesh->y, 0.0, mesh->cubature->op_tmp[2]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB_VDS, mesh->y, 0.0, mesh->cubature->op_tmp[3]);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB_V, factor, 0.0, cFactor);

  op_par_loop(poisson_cubature_op, "poisson_cubature_op", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::CUB_VDR), DG_ORDER * DG_CUB_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::CUB_VDS), DG_ORDER * DG_CUB_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[0], -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[1], -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[2], -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->op_tmp[3], -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(mesh->cubature->J, -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(cFactor, -1, OP_ID, DG_CUB_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_WRITE));
}

void PoissonMatrix2D::calc_gauss_sub_mat() {
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, factor, 0.0, gFactor);

  // TODO change for edges with different orders on each side
  op_par_loop(poisson_gauss_grad, "poisson_gauss_grad", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DR), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DS), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DR), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DS), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DR), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DS), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP0), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP1), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP2), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->x, -2, mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -2, mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(h,       -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(gFactor, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(op1, 0, mesh->face2cells, DG_NP * DG_NP, "double", OP_INC),
              op_arg_dat(op1, 1, mesh->face2cells, DG_NP * DG_NP, "double", OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_WRITE));

  // If not dirichlet BC, kernel will assume it is a neumann bc
  if(mesh->bface2cells) {
    op_par_loop(poisson_gauss_grad_b, "poisson_gauss_grad_b", mesh->bfaces,
                op_arg_dat(mesh->order, 0, mesh->bface2cells, 1, "int", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DR), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F0DS), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DR), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F1DS), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DR), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_F2DS), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP0), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP1), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_gbl(constants->get_mat_ptr(DGConstants::GAUSS_FINTERP2), DG_ORDER * DG_GF_NP * DG_NP, "double", OP_READ),
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
                op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
                op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, "double", OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->sJ, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(h,       0, mesh->bface2cells, 1, "double", OP_READ),
                op_arg_dat(gFactor, 0, mesh->bface2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(op1, 0, mesh->bface2cells, DG_NP * DG_NP, "double", OP_INC),
                op_arg_dat(op_bc, -1, OP_ID, DG_GF_NP * DG_NP, "double", OP_WRITE));
  }
}

void PoissonMatrix2D::calc_mm_mat() {
  op_par_loop(poisson_mm, "poisson_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mmFactor, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_INC));
}
