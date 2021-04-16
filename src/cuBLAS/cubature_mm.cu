#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_cubature_mm(cublasHandle_t handle, const int numCells,
                               const double *temp_d, double *mm_d) {
  // double *cubV_d;
  // cudaMalloc((void**)&cubV_d, 46 * 15 * sizeof(double));
  // cudaMemcpy(cubV_d, cubV, 46 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  for(int i = 0; i < numCells; i++) {
    const double *temp = temp_d + i * 46 * 15;
    double *mm = mm_d + i * 15 * 15;

    double alpha = 1.0;
    double beta = 0.0;
    cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, 15, 15, 46, &alpha, constants->cubV_d, 15, temp, 15, &beta, mm, 15);
  }

  // cudaFree(cubV_d);
}

void cubature_mm_blas(INSData *nsData, CubatureData *cubData) {
  // Initialise cuBLAS
  // cublasHandle_t handle;
  // cublasCreate(&handle);
  // cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg mm_cubature_args[] = {
    op_arg_dat(cubData->temp, -1, OP_ID, 46 * 15, "double", OP_READ),
    op_arg_dat(cubData->mm, -1, OP_ID, 15 * 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 2, mm_cubature_args);

  cublas_cubature_mm(constants->handle, nsData->numCells, (double *)cubData->temp->data_d,
                     (double *)cubData->mm->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(2, mm_cubature_args);
  // Free resources used by cuBLAS
  // cublasDestroy(handle);
}
