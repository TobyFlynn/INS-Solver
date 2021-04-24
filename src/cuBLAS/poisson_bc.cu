#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_poisson_bc(cublasHandle_t handle, const int numCells,
                              const double *temp0_d, const double *temp1_d,
                              double *gx_d, double *gy_d) {
  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double *mat_d;
  cudaMalloc((void**)&mat_d, 21 * 15 * sizeof(double));
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, 21, 15, &alpha, constants->invMass_d, 15, constants->gInterp_d, 15, &beta, mat_d, 15);

  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 21, &alpha, mat_d, 15, temp0_d, 21, &beta, gx_d, 15);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 21, &alpha, mat_d, 15, temp1_d, 21, &beta, gy_d, 15);

  cudaFree(mat_d);
}

void poisson_bc_blas(INSData *data, Poisson_MF *poisson) {
  // Make sure OP2 data is in the right place
  op_arg poisson_args[] = {
    op_arg_dat(poisson->uFluxX, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->uFluxY, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->gradx, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(poisson->grady, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 4, poisson_args);

  cublas_poisson_bc(constants->handle, data->numCells, (double *)poisson->uFluxX->data_d,
                    (double *)poisson->uFluxY->data_d, (double *)poisson->gradx->data_d,
                    (double *)poisson->grady->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(4, poisson_args);
}
