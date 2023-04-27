#include "matrices/3d/poisson_coarse_matrix_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

PoissonCoarseMatrix3D::PoissonCoarseMatrix3D(DGMesh3D *m) {
  mesh = m;
  _mesh = m;

  petscMatInit = false;

  op1    = op_decl_dat(mesh->cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_op1");
  op2[0] = op_decl_dat(mesh->faces, DG_NP_N1 * DG_NP_N1, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_op20");
  op2[1] = op_decl_dat(mesh->faces, DG_NP_N1 * DG_NP_N1, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_op21");
  opbc   = op_decl_dat(mesh->bfaces, DG_NP_N1 * DG_NPF_N1, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_opbc");

  glb_ind = op_decl_dat(mesh->cells, 1, "int", (int *)NULL, "poisson_matrix_glb_ind");
  glb_indL = op_decl_dat(mesh->faces, 1, "int", (int *)NULL, "poisson_matrix_glb_indL");
  glb_indR = op_decl_dat(mesh->faces, 1, "int", (int *)NULL, "poisson_matrix_glb_indR");
  orderL  = op_decl_dat(mesh->faces, 1, "int", (int *)NULL, "poisson_orderL");
  orderR  = op_decl_dat(mesh->faces, 1, "int", (int *)NULL, "poisson_orderR");
}

void PoissonCoarseMatrix3D::calc_mat() {
  timer->startTimer("PoissonCoarseMatrix3D - calc_mat");
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  petscMatResetRequired = true;
  timer->endTimer("PoissonCoarseMatrix3D - calc_mat");
}

void PoissonCoarseMatrix3D::apply_bc(op_dat rhs, op_dat bc) {
  timer->startTimer("PoissonCoarseMatrix3D - apply_bc");
  if(mesh->bface2cells) {
    op_par_loop(poisson_coarse_matrix_3d_apply_bc, "poisson_coarse_matrix_3d_apply_bc", mesh->bfaces,
                op_arg_dat(opbc, -1, OP_ID, DG_NPF_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
                op_arg_dat(bc,   -1, OP_ID, DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(rhs,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
  timer->endTimer("PoissonCoarseMatrix3D - apply_bc");
}

void PoissonCoarseMatrix3D::calc_op1() {
  timer->startTimer("PoissonCoarseMatrix3D - calc_op1");
  op_par_loop(poisson_coarse_matrix_3d_op1, "poisson_coarse_matrix_3d_op1", mesh->cells,
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
              op_arg_dat(op1, -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonCoarseMatrix3D - calc_op1");
}

void PoissonCoarseMatrix3D::calc_op2() {
  timer->startTimer("PoissonCoarseMatrix3D - calc_op2");
  // TODO full p-adaptivity
  op_par_loop(poisson_coarse_matrix_3d_op2, "poisson_coarse_matrix_3d_op2", mesh->faces,
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
              op_arg_dat(op1, 0, mesh->face2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC),
              op_arg_dat(op1, 1, mesh->face2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonCoarseMatrix3D - calc_op2");
}

void PoissonCoarseMatrix3D::calc_opbc() {
  timer->startTimer("PoissonCoarseMatrix3D - calc_opbc");
  if(mesh->bface2cells) {
    op_par_loop(poisson_coarse_matrix_3d_bop, "poisson_coarse_matrix_3d_bop", mesh->bfaces,
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
                op_arg_dat(op1, 0, mesh->bface2cells, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_INC));

    op_par_loop(poisson_coarse_matrix_3d_opbc, "poisson_coarse_matrix_3d_opbc", mesh->bfaces,
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
                op_arg_dat(opbc, -1, OP_ID, DG_NPF_N1 * DG_NP_N1, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("PoissonCoarseMatrix3D - calc_opbc");
}

void PoissonCoarseMatrix3D::calc_glb_ind() {
  timer->startTimer("PoissonCoarseMatrix3D - calc_glb_ind");
  set_glb_ind();
  op_par_loop(copy_to_edges, "copy_to_edges", mesh->faces,
              op_arg_dat(glb_ind, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));

  op_par_loop(copy_to_edges, "copy_to_edges", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_WRITE));
  timer->endTimer("PoissonCoarseMatrix3D - calc_glb_ind");
}
