#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

void op2_gemv_batch(bool transpose, int m, int n, double alpha, op_dat a, int lda, op_dat x, double beta, op_dat y) {
  op_arg gemv_args[] = {
    op_arg_dat(a, -1, OP_ID, a->dim, "double", OP_READ),
    op_arg_dat(x, -1, OP_ID, x->dim, "double", OP_READ),
    op_arg_dat(y, -1, OP_ID, y->dim, "double", OP_RW)
  };
  op_mpi_halo_exchanges(a->set, 3, gemv_args);

  int setSize = a->set->size;

  if(transpose) {
    for(int i = 0; i < setSize; i++) {
      const double *a_c = (double *)a->data + i * a->dim;
      const double *x_c = (double *)x->data + i * x->dim;
      double *y_c = (double *)y->data + i * y->dim;

      cblas_dgemv(CblasColMajor, CblasTrans, m, n, alpha, a_c, lda, x_c, 1, beta, y_c, 1);
    }
  } else {
    for(int i = 0; i < setSize; i++) {
      const double *a_c = (double *)a->data + i * a->dim;
      const double *x_c = (double *)x->data + i * x->dim;
      double *y_c = (double *)y->data + i * y->dim;

      cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, alpha, a_c, lda, x_c, 1, beta, y_c, 1);
    }
  }

  op_mpi_set_dirtybit(3, gemv_args);
}
