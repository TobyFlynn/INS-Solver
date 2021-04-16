#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_gauss_interp(cublasHandle_t handle, const int numCells,
                                const double *in_d, double *out_d) {
  // double *interp_d;
  // cudaMalloc((void**)&interp_d, 21 * 15 * sizeof(double));
  // cudaMemcpy(interp_d, gInterp, 21 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 21, numCells, 15, &alpha, constants->gInterp_d, 15, in_d, 15, &beta, out_d, 21);
  // cudaFree(interp_d);
}

void gauss_interp_blas(INSData *data, op_dat input, op_dat output) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg interp_args[] = {
    op_arg_dat(input, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(output, -1, OP_ID, 21, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 2, interp_args);

  cublas_gauss_interp(handle, data->numCells, (double *)input->data_d,
                      (double *)output->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(2, interp_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
