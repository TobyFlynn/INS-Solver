#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_init_grid(cublasHandle_t handle, const int numCells,
                        const double *node_coords, const int *cell2nodes,
                        double *x_d, double *y_d, double *xr_d, double *xs_d,
                        double *yr_d, double *ys_d) {
  double *ones_d;
  cudaMalloc((void**)&ones_d, 15 * sizeof(double));
  cudaMemcpy(ones_d, ones, 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *r_d;
  cudaMalloc((void**)&r_d, 15 * sizeof(double));
  cudaMemcpy(r_d, r, 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *s_d;
  cudaMalloc((void**)&s_d, 15 * sizeof(double));
  cudaMemcpy(s_d, s, 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *temp_d;
  cudaMalloc((void**)&temp_d, numCells * 15 * sizeof(double));

  double *Dr_d;
  cudaMalloc((void**)&Dr_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Dr_d, Dr, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *Ds_d;
  cudaMalloc((void**)&Ds_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Ds_d, Ds, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  for(int c = 0; c < numCells; c++) {
    // Get nodes for this cell (on host)
    const double *n0 = &node_coords[2 * cell2nodes[3 * c]];
    const double *n1 = &node_coords[2 * cell2nodes[3 * c + 1]];
    const double *n2 = &node_coords[2 * cell2nodes[3 * c + 2]];

    double *temp = temp_d + c * 15;
    double *x = x_d + c * 15;
    double *y = y_d + c * 15;

    double alpha = 1.0;
    cublasDcopy(handle, 15, ones_d, 1, x, 1);
    cublasDaxpy(handle, 15, &alpha, r_d, 1, x, 1);
    alpha = 0.5 * n1[0];
    cublasDscal(handle, 15, &alpha, x, 1);
    cublasDcopy(handle, 15, ones_d, 1, temp, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, s_d, 1, temp, 1);
    alpha = 0.5 * n2[0];
    cublasDaxpy(handle, 15, &alpha, temp, 1, x, 1);
    cublasDcopy(handle, 15, s_d, 1, temp, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, r_d, 1, temp, 1);
    alpha = -0.5 * n0[0];
    cublasDaxpy(handle, 15, &alpha, temp, 1, x, 1);

    cublasDcopy(handle, 15, ones_d, 1, y, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, r_d, 1, y, 1);
    alpha = 0.5 * n1[1];
    cublasDscal(handle, 15, &alpha, y, 1);
    cublasDcopy(handle, 15, ones_d, 1, temp, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, s_d, 1, temp, 1);
    alpha = 0.5 * n2[1];
    cublasDaxpy(handle, 15, &alpha, temp, 1, y, 1);
    cublasDcopy(handle, 15, s_d, 1, temp, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, r_d, 1, temp, 1);
    alpha = -0.5 * n0[1];
    cublasDaxpy(handle, 15, &alpha, temp, 1, y, 1);
  }

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha2 = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, Dr_d, 15, x_d, 15, &beta, xr_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, Ds_d, 15, x_d, 15, &beta, xs_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, Dr_d, 15, y_d, 15, &beta, yr_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, Ds_d, 15, y_d, 15, &beta, ys_d, 15);

  cudaFree(ones_d);
  cudaFree(r_d);
  cudaFree(temp_d);
  cudaFree(Dr_d);
  cudaFree(Ds_d);
}

void init_grid_blas(INSData *nsData) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg init_grid_args[] = {
    op_arg_dat(nsData->x, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->y, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->xr, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->xs, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->yr, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->ys, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 6, init_grid_args);

  cublas_init_grid(handle, nsData->numCells, (double *)nsData->node_coords->data,
                   (int *)nsData->cell2nodes->map, (double *)nsData->x->data_d,
                   (double *)nsData->y->data_d, (double *)nsData->xr->data_d,
                   (double *)nsData->xs->data_d, (double *)nsData->yr->data_d,
                   (double *)nsData->ys->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(6, init_grid_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
