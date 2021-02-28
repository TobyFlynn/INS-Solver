#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_pressure_rhs(cublasHandle_t handle, const int numCells,
                                const double *div_d, const double *dPdN_d,
                                double *rhs_d) {
  double *LIFT_d;
  cudaMalloc((void**)&LIFT_d, 15 * 15 * sizeof(double));
  cudaMemcpy(LIFT_d, LIFT, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *MASS_d;
  cudaMalloc((void**)&MASS_d, 15 * 15 * sizeof(double));
  cudaMemcpy(MASS_d, MASS, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta1 = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, MASS_d, 15, div_d, 15, &beta1, rhs_d, 15);
  double alpha2 = -1.0;
  double beta2 = -1.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, LIFT_d, 15, dPdN_d, 15, &beta2, rhs_d, 15);

  cudaFree(LIFT_d);
  cudaFree(MASS_d);
}

void pressure_rhs_blas(INSData *nsData, int ind) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg pressure_rhs_args[] = {
    op_arg_dat(nsData->divVelT, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->dPdN[(ind + 1) % 2], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->pRHS, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 3, pressure_rhs_args);

  cublas_pressure_rhs(handle, nsData->numCells, (double *)nsData->divVelT->data_d,
                      (double *)nsData->dPdN[(ind + 1) % 2]->data_d,
                      (double *)nsData->pRHS->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(3, pressure_rhs_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
