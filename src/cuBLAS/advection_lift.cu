#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_advection_lift(cublasHandle_t handle, const int numCells,
                                  const double *flux0_d, const double *flux1_d,
                                  double *N0_d, double *N1_d) {
  // double *LIFT_d;
  // cudaMalloc((void**)&LIFT_d, 15 * 15 * sizeof(double));
  // cudaMemcpy(LIFT_d, LIFT, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 1.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, constants->LIFT_d, 15, flux0_d, 15, &beta, N0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, constants->LIFT_d, 15, flux1_d, 15, &beta, N1_d, 15);
  // cudaFree(LIFT_d);
}

void advection_lift_blas(INSData *nsData, int ind) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg advec_lift_args[] = {
    op_arg_dat(nsData->flux[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->flux[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->N[ind][0], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(nsData->N[ind][1], -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 4, advec_lift_args);

  cublas_advection_lift(handle, nsData->numCells, (double *)nsData->flux[0]->data_d,
                     (double *)nsData->flux[1]->data_d, (double *)nsData->N[ind][0]->data_d,
                     (double *)nsData->N[ind][1]->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(4, advec_lift_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
