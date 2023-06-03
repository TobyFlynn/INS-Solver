#include "matrices/2d/poisson_coarse_matrix_2d.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"

extern Timing *timer;
extern DGConstants *constants;

PoissonCoarseMatrix2D::PoissonCoarseMatrix2D(DGMesh2D *m) {
  mesh = m;
  _mesh = m;
  petscMatInit = false;

  op1    = op_decl_dat(mesh->cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, (DG_FP *)NULL, "poisson_op1");
  op2[0] = op_decl_dat(mesh->faces, DG_NP_N1 * DG_NP_N1, DG_FP_STR, (DG_FP *)NULL, "poisson_op20");
  op2[1] = op_decl_dat(mesh->faces, DG_NP_N1 * DG_NP_N1, DG_FP_STR, (DG_FP *)NULL, "poisson_op21");
  opbc   = op_decl_dat(mesh->bfaces, DG_NUM_FACES * DG_NPF_N1 * DG_NP_N1, DG_FP_STR, (DG_FP *)NULL, "poisson_opbc");

  glb_ind   = op_decl_dat(mesh->cells, 1, "int", (int *)NULL, "poisson_glb_ind");
  glb_indL  = op_decl_dat(mesh->faces, 1, "int", (int *)NULL, "poisson_glb_indL");
  glb_indR  = op_decl_dat(mesh->faces, 1, "int", (int *)NULL, "poisson_glb_indR");
}

PoissonCoarseMatrix2D::~PoissonCoarseMatrix2D() {
  if(petscMatInit)
    MatDestroy(&pMat);
}

void PoissonCoarseMatrix2D::calc_mat() {
  timer->startTimer("PoissonCoarseMatrix2D - calc_mat");
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  petscMatResetRequired = true;
  timer->endTimer("PoissonCoarseMatrix2D - calc_mat");
}

void PoissonCoarseMatrix2D::calc_glb_ind() {
  timer->startTimer("PoissonCoarseMatrix2D - calc_glb_ind");
  set_glb_ind();
  op_par_loop(copy_to_edges, "copy_to_edges", mesh->faces,
              op_arg_dat(glb_ind, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));
  timer->endTimer("PoissonCoarseMatrix2D - calc_glb_ind");
}

void PoissonCoarseMatrix2D::apply_bc(op_dat rhs, op_dat bc) {
  timer->startTimer("PoissonCoarseMatrix2D - apply_bc");
  /*
  if(mesh->bface2cells) {
    op_par_loop(poisson_coarse_matrix_2d_apply_bc, "poisson_coarse_matrix_2d_apply_bc", mesh->bfaces,
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(opbc, -1, OP_ID, DG_GF_NP * DG_NP_N1, DG_FP_STR, OP_READ),
                op_arg_dat(bc,    0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rhs,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
  */
  timer->endTimer("PoissonCoarseMatrix2D - apply_bc");
}

void PoissonCoarseMatrix2D::calc_op1() {
  timer->startTimer("PoissonCoarseMatrix2D - calc_op1");
  op_par_loop(poisson_coarse_matrix_2d_op1, "poisson_coarse_matrix_2d_op1", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonCoarseMatrix2D - calc_op1");
}

void PoissonCoarseMatrix2D::calc_op2() {
  timer->startTimer("PoissonCoarseMatrix2D - calc_op2");
  // TODO write this matrix (need new constant matrices from toolkit)
  /*
  op_par_loop(poisson_coarse_matrix_2d_op2, "poisson_coarse_matrix_2d_op2", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -2, mesh->face2cells, 5, DG_FP_STR, OP_READ),
              op_arg_dat(op1, 0, mesh->face2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC),
              op_arg_dat(op1, 1, mesh->face2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE));
  */
  timer->endTimer("PoissonCoarseMatrix2D - calc_op2");
}

void PoissonCoarseMatrix2D::calc_opbc() {
  timer->startTimer("PoissonCoarseMatrix2D - calc_opbc");
  /*
  if(mesh->bface2cells) {
    op_par_loop(poisson_coarse_matrix_2d_bop, "poisson_coarse_matrix_2d_bop", mesh->bfaces,
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
                op_arg_dat(op1, 0, mesh->bface2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC),
                op_arg_dat(opbc, -1, OP_ID, DG_GF_NP * DG_NP_N1, DG_FP_STR, OP_WRITE));
  }
  */
  timer->endTimer("PoissonCoarseMatrix2D - calc_opbc");
}
