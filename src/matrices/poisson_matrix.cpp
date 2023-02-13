#include "matrices/poisson_matrix.h"

#include "op_seq.h"

#include "timing.h"

extern Timing *timer;

void PoissonMatrix::set_bc_types(op_dat bc_ty) {
  bc_types = bc_ty;
}

bool PoissonMatrix::getPETScMat(Mat** mat) {
  bool reset = false;
  if(petscMatResetRequired) {
    setPETScMatrix();
    petscMatResetRequired = false;
    *mat = &pMat;
    reset = true;
  }

  return reset;
}

void PoissonMatrix::mult(op_dat in, op_dat out) {
  timer->startTimer("M - Mult");
  timer->startTimer("M - Mult - Cells");
  op_par_loop(poisson_mult_cells, "poisson_mult_cells", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_WRITE));
  timer->endTimer("M - Mult - Cells");
  timer->startTimer("M - Mult - Faces");
  op_par_loop(poisson_mult_faces, "poisson_mult_faces", _mesh->faces,
              op_arg_dat(_mesh->order, 0, _mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(in,      0, _mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out,     0, _mesh->face2cells, DG_NP, "double", OP_INC),
              op_arg_dat(_mesh->order, 1, _mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(in,      1, _mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out,     1, _mesh->face2cells, DG_NP, "double", OP_INC));
  timer->endTimer("M - Mult - Faces");
  timer->endTimer("M - Mult");
}

void PoissonMatrix::multJacobi(op_dat in, op_dat out) {
  mult(in, out);

  op_par_loop(poisson_mult_jacobi, "poisson_mult_jacobi", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_RW));
}
