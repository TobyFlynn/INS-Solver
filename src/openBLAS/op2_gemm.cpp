#include "cblas.h"

#include "op_seq.h"
#include "blas_calls.h"

void op2_gemm(bool transposeA, bool transposeB, int m, int n, int k, double alpha, double *a_ptr, int lda, op_dat b, int ldb, double beta, op_dat c, int ldc) {
  op_arg gemv_args[] = {
    op_arg_dat(b, -1, OP_ID, b->dim, "double", OP_READ),
    op_arg_dat(c, -1, OP_ID, c->dim, "double", OP_RW)
  };
  op_mpi_halo_exchanges(b->set, 2, gemv_args);

  int setSize = b->set->size;

  if(transposeA) {
    if(transposeB) {
      for(int i = 0; i < setSize; i++) {
        const double *b_c = (double *)b->data + i * b->dim;
        double *c_c = (double *)c->data + i * c->dim;

        cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, m, n, k, alpha, a_ptr, lda, b_c, ldb, beta, c_c, ldc);
      }
    } else {
      for(int i = 0; i < setSize; i++) {
        const double *b_c = (double *)b->data + i * b->dim;
        double *c_c = (double *)c->data + i * c->dim;

        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, a_ptr, lda, b_c, ldb, beta, c_c, ldc);
      }
    }
  } else {
    if(transposeB) {
      for(int i = 0; i < setSize; i++) {
        const double *b_c = (double *)b->data + i * b->dim;
        double *c_c = (double *)c->data + i * c->dim;

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, a_ptr, lda, b_c, ldb, beta, c_c, ldc);
      }
    } else {
      for(int i = 0; i < setSize; i++) {
        const double *b_c = (double *)b->data + i * b->dim;
        double *c_c = (double *)c->data + i * c->dim;

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a_ptr, lda, b_c, ldb, beta, c_c, ldc);
      }
    }
  }

  op_mpi_set_dirtybit(2, gemv_args);
}

void op2_gemm(bool transposeA, bool transposeB, int m, int n, int k, double alpha, op_dat a, int lda, double *b_ptr, int ldb, double beta, op_dat c, int ldc) {
  op_arg gemv_args[] = {
    op_arg_dat(a, -1, OP_ID, a->dim, "double", OP_READ),
    op_arg_dat(c, -1, OP_ID, c->dim, "double", OP_RW)
  };
  op_mpi_halo_exchanges(a->set, 2, gemv_args);

  int setSize = a->set->size;

  if(transposeA) {
    if(transposeB) {
      for(int i = 0; i < setSize; i++) {
        const double *a_c = (double *)a->data + i * a->dim;
        double *c_c = (double *)c->data + i * c->dim;

        cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, m, n, k, alpha, a_c, lda, b_ptr, ldb, beta, c_c, ldc);
      }
    } else {
      for(int i = 0; i < setSize; i++) {
        const double *a_c = (double *)a->data + i * a->dim;
        double *c_c = (double *)c->data + i * c->dim;

        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, a_c, lda, b_ptr, ldb, beta, c_c, ldc);
      }
    }
  } else {
    if(transposeB) {
      for(int i = 0; i < setSize; i++) {
        const double *a_c = (double *)a->data + i * a->dim;
        double *c_c = (double *)c->data + i * c->dim;

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, a_c, lda, b_ptr, ldb, beta, c_c, ldc);
      }
    } else {
      for(int i = 0; i < setSize; i++) {
        const double *a_c = (double *)a->data + i * a->dim;
        double *c_c = (double *)c->data + i * c->dim;

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a_c, lda, b_ptr, ldb, beta, c_c, ldc);
      }
    }
  }

  op_mpi_set_dirtybit(2, gemv_args);
}
