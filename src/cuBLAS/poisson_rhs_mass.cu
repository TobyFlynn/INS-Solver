#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_poisson_rhs_mass(cublasHandle_t handle, const int numCells, const double *mm_d,
                                    const double *u_d, double *rhs_d, double factor) {
  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = factor;
  double beta = 1.0;

  for(int i = 0; i < numCells; i++) {
    const double *mm = mm_d + i * 15 * 15;
    const double *u = u_d + i * 15;
    double *rhs = rhs_d + i * 15;

    cublasDgemv(handle, CUBLAS_OP_N, 15, 15, &alpha, mm, 15, u, 1, &beta, rhs, 1);
  }
}

void poisson_rhs_mass_blas(INSData *data, CubatureData *cubatureData, Poisson_MF *poisson, double factor) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg mass_args[] = {
    op_arg_dat(cubatureData->mm, -1, OP_ID, 15 * 15, "double", OP_READ),
    op_arg_dat(poisson->u, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(poisson->rhs, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 3, mass_args);

  cublas_poisson_rhs_mass(handle, data->numCells, (double *)cubatureData->mm->data_d,
                          (double *)poisson->u->data_d, (double *)poisson->rhs->data_d, factor);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(3, mass_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
