//
// auto-generated by op2.py
//

#include "cublas_v2.h"

#include  "op_lib_cpp.h"

//
// op_par_loop declarations
//
#ifdef OPENACC
#ifdef __cplusplus
extern "C" {
#endif
#endif
#ifdef OPENACC
#ifdef __cplusplus
}
#endif
#endif

#include "../blas_calls.h"

inline void cublas_viscosity_rhs(cublasHandle_t handle, const int numCells,
                                 const double *qtt0_d, const double *qtt1_d,
                                 double *vis0_d, double *vis1_d) {
  // double *visMat_d;
  // cudaMalloc((void**)&visMat_d, 15 * 15 * sizeof(double));
  // cudaMemcpy(visMat_d, visMat, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *MASS_d;
  cudaMalloc((void**)&MASS_d, 15 * 15 * sizeof(double));
  cudaMemcpy(MASS_d, MASS, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  // cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, visMat_d, 15, qtt0_d, 15, &beta, temp0_d, 15);
  // cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, visMat_d, 15, qtt1_d, 15, &beta, temp1_d, 15);

  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, MASS_d, 15, qtt0_d, 15, &beta, vis0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, MASS_d, 15, qtt1_d, 15, &beta, vis1_d, 15);

  // cudaFree(visMat_d);
  cudaFree(MASS_d);
}

void viscosity_rhs_blas(INSData *nsData) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg viscosity_rhs_args[] = {
    op_arg_dat(nsData->QTT[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->QTT[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->visRHS[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->visRHS[1], -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 4, viscosity_rhs_args);

  cublas_viscosity_rhs(handle, nsData->numCells, (double *)nsData->QTT[0]->data_d,
                      (double *)nsData->QTT[1]->data_d, (double *)nsData->visRHS[0]->data_d,
                      (double *)nsData->visRHS[1]->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(4, viscosity_rhs_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}