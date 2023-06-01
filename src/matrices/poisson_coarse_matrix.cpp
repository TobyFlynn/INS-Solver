#include "matrices/poisson_coarse_matrix.h"

#include "op_seq.h"

#include "timing.h"

extern Timing *timer;

void PoissonCoarseMatrix::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonCoarseMatrix - mult");
  timer->startTimer("PoissonCoarseMatrix - mult Cells");
  op_par_loop(poisson_mult_cells_coarse, "poisson_mult_cells_coarse", _mesh->cells,
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonCoarseMatrix - mult Cells");
  timer->startTimer("PoissonCoarseMatrix - mult Faces");
  op_par_loop(poisson_mult_faces_coarse, "poisson_mult_faces_coarse", _mesh->faces,
              op_arg_dat(in,      0, _mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
              op_arg_dat(out,     0, _mesh->face2cells, DG_NP, DG_FP_STR, OP_INC),
              op_arg_dat(in,      1, _mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
              op_arg_dat(out,     1, _mesh->face2cells, DG_NP, DG_FP_STR, OP_INC));
  timer->endTimer("PoissonCoarseMatrix - mult Faces");
  timer->endTimer("PoissonCoarseMatrix - mult");
}

void PoissonCoarseMatrix::multJacobi(op_dat in, op_dat out) {
  timer->startTimer("PoissonCoarseMatrix - multJacobi");
  mult(in, out);

  op_par_loop(poisson_mult_jacobi_coarse, "poisson_mult_jacobi_coarse", _mesh->cells,
              op_arg_dat(op1, -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("PoissonCoarseMatrix - multJacobi");
}

#ifdef INS_BUILD_WITH_AMGX
bool PoissonCoarseMatrix::getAmgXMat(AMGX_matrix_handle** mat) {
  timer->startTimer("PoissonCoarseMatrix - getAmgXMat");
  bool reset = false;
  if(petscMatResetRequired) {
    setAmgXMatrix();
    petscMatResetRequired = false;
    *mat = &amgx_mat;
    reset = true;
  }
  timer->endTimer("PoissonCoarseMatrix - getAmgXMat");

  return reset;
}
#endif

#ifdef INS_BUILD_WITH_HYPRE
bool PoissonCoarseMatrix::getHYPREMat(HYPRE_ParCSRMatrix** mat) {
  timer->startTimer("PoissonCoarseMatrix - getHYPREMat");
  bool reset = false;
  if(petscMatResetRequired) {
    setHYPREMatrix();
    HYPRE_IJMatrixGetObject(hypre_mat, (void**)&hypre_parcsr_mat);
    petscMatResetRequired = false;
    reset = true;
  }
  *mat = &hypre_parcsr_mat;
  timer->endTimer("PoissonCoarseMatrix - getHYPREMat");

  return reset;
}

void PoissonCoarseMatrix::getHYPRERanges(int *ilower, int *iupper, int *jlower, int *jupper) {
  HYPRE_IJMatrixGetLocalRange(hypre_mat, ilower, iupper, jlower, jupper);
}
#endif
