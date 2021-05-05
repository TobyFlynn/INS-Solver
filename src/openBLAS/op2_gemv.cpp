#include "cblas.h"

#include "op_seq.h"
#include "blas_calls.h"

void op2_gemv(bool transpose, int m, int n, double alpha, double *A_ptr, int lda, op_dat x, double beta, op_dat y) {
  op_arg gemv_args[] = {
    op_arg_dat(x, -1, OP_ID, x->dim, "double", OP_READ),
    op_arg_dat(y, -1, OP_ID, y->dim, "double", OP_RW)
  };
  op_mpi_halo_exchanges(x->set, 2, gemv_args);

  int setSize = x->set->size;

  if(transpose) {
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, setSize, n, alpha, A_ptr, lda, (double *)x->data, x->dim, beta, (double *)y->data, y->dim);
  } else {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, setSize, n, alpha, A_ptr, lda, (double *)x->data, x->dim, beta, (double *)y->data, y->dim);
  }

  op_mpi_set_dirtybit(2, gemv_args);
}
