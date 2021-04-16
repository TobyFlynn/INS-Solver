#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_gauss_op(cublasHandle_t handle, const int numCells,
                            double *op_d, const double *mD_d, const double *term0_d,
                            const double *term1_d, const double *term2_d, int face) {
  double *gFInterp_d;
  // cudaMalloc((void**)&gFInterp_d, 7 * 15 * sizeof(double));
  if(face == 0) {
    // cudaMemcpy(gFInterp_d, gFInterp0, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
    gFInterp_d = constants->gFInterp0_d;
  } else if (face == 1) {
    // cudaMemcpy(gFInterp_d, gFInterp1, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
    gFInterp_d = constants->gFInterp1_d;
  } else {
    // cudaMemcpy(gFInterp_d, gFInterp2, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
    gFInterp_d = constants->gFInterp2_d;
  }

  for(int c = 0; c < numCells; c++) {
    double *op = op_d + c * 15 * 15;
    const double *mD = mD_d + c * 7 * 15;
    const double *term0 = term0_d + c * 7 * 15;
    const double *term1 = term1_d + c * 7 * 15;
    const double *term2 = term2_d + c * 7 * 15;

    double alpha = 1.0;
    double beta = 0.0;
    cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_T, 15, 15, 7, &alpha, term0, 7, gFInterp_d, 15, &beta, op, 15);
    double alpha2 = -1.0;
    double beta2 = 1.0;
    cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_T, 15, 15, 7, &alpha2, term1, 7, mD, 15, &beta2, op, 15);
    cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_T, 15, 15, 7, &alpha2, term2, 7, gFInterp_d, 15, &beta2, op, 15);
  }

  // cudaFree(gFInterp_d);
}

void gauss_op_blas(INSData *nsData, GaussData *gaussData) {
  // Initialise cuBLAS
  // cublasHandle_t handle;
  // cublasCreate(&handle);
  // cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg gauss_args[] = {
    // Face 0
    op_arg_dat(gaussData->OP[0], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->mD[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    // Face 1
    op_arg_dat(gaussData->OP[1], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->mD[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    // Face 2
    op_arg_dat(gaussData->OP[2], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->mD[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[2], -1, OP_ID, 7 * 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 15, gauss_args);

  cublas_gauss_op(constants->handle, nsData->numCells, (double *)gaussData->OP[0]->data_d,
                   (double *)gaussData->mD[0]->data_d, (double *)gaussData->mDx[0]->data_d,
                   (double *)gaussData->mDx[1]->data_d, (double *)gaussData->mDx[2]->data_d, 0);

  cublas_gauss_op(constants->handle, nsData->numCells, (double *)gaussData->OP[1]->data_d,
                   (double *)gaussData->mD[1]->data_d, (double *)gaussData->mDy[0]->data_d,
                   (double *)gaussData->mDy[1]->data_d, (double *)gaussData->mDy[2]->data_d, 1);

  cublas_gauss_op(constants->handle, nsData->numCells, (double *)gaussData->OP[2]->data_d,
                   (double *)gaussData->mD[2]->data_d, (double *)gaussData->pDx[0]->data_d,
                   (double *)gaussData->pDx[1]->data_d, (double *)gaussData->pDx[2]->data_d, 2);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(15, gauss_args);
  // Free resources used by cuBLAS
  // cublasDestroy(handle);
}
