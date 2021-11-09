#include "cublas_v2.h"

#include "op_seq.h"
#include "blas_calls.h"

inline void cublas_init_gauss_grad_neighbour(cublasHandle_t handle, const int numCells, const int *reverse,
                        const double *x_d, const double *y_d, double *gxr_d,
                        double *gxs_d, double *gyr_d, double *gys_d) {
  // Calc Grad Matrices
  double alpha = 1.0;
  double beta = 0.0;
  for(int c = 0; c < numCells; c++) {
    const double *x = x_d + c * 10;
    const double *y = y_d + c * 10;
    double *gxr = gxr_d + c * 18;
    double *gxs = gxs_d + c * 18;
    double *gyr = gyr_d + c * 18;
    double *gys = gys_d + c * 18;

    // Face 0
    if(reverse[3 * c]) {
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0DrR_d, 6, x, 1, &beta, gxr, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0DsR_d, 6, x, 1, &beta, gxs, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0DrR_d, 6, y, 1, &beta, gyr, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0DsR_d, 6, y, 1, &beta, gys, 1);
    } else {
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0Dr_d, 6, x, 1, &beta, gxr, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0Ds_d, 6, x, 1, &beta, gxs, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0Dr_d, 6, y, 1, &beta, gyr, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0Ds_d, 6, y, 1, &beta, gys, 1);
    }

    // Face 1
    if(reverse[3 * c + 1]) {
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1DrR_d, 6, x, 1, &beta, gxr + 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1DsR_d, 6, x, 1, &beta, gxs + 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1DrR_d, 6, y, 1, &beta, gyr + 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1DsR_d, 6, y, 1, &beta, gys + 6, 1);
    } else {
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1Dr_d, 6, x, 1, &beta, gxr + 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1Ds_d, 6, x, 1, &beta, gxs + 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1Dr_d, 6, y, 1, &beta, gyr + 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1Ds_d, 6, y, 1, &beta, gys + 6, 1);
    }

    // Face 2
    if(reverse[3 * c + 2]) {
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2DrR_d, 6, x, 1, &beta, gxr + 2 * 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2DsR_d, 6, x, 1, &beta, gxs + 2 * 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2DrR_d, 6, y, 1, &beta, gyr + 2 * 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2DsR_d, 6, y, 1, &beta, gys + 2 * 6, 1);
    } else {
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2Dr_d, 6, x, 1, &beta, gxr + 2 * 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2Ds_d, 6, x, 1, &beta, gxs + 2 * 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2Dr_d, 6, y, 1, &beta, gyr + 2 * 6, 1);
      cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2Ds_d, 6, y, 1, &beta, gys + 2 * 6, 1);
    }
  }
}

void init_gauss_grad_neighbour_blas(DGMesh *mesh, INSData *data) {
  // Make sure OP2 data is in the right place
  op_arg init_grad_args[] = {
    op_arg_dat(mesh->x, -1, OP_ID, 10, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, 10, "double", OP_READ),
    op_arg_dat(data->reverse, -1, OP_ID, 3, "int", OP_READ),
    op_arg_dat(data->grx, -1, OP_ID, 18, "double", OP_WRITE),
    op_arg_dat(data->gsx, -1, OP_ID, 18, "double", OP_WRITE),
    op_arg_dat(data->gry, -1, OP_ID, 18, "double", OP_WRITE),
    op_arg_dat(data->gsy, -1, OP_ID, 18, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 7, init_grad_args);

  int setSize = mesh->x->set->size;

  int *reverse = (int *)malloc(3 * setSize * sizeof(int));
  cudaMemcpy(reverse, data->reverse->data_d, 3 * setSize * sizeof(int), cudaMemcpyDeviceToHost);

  cublas_init_gauss_grad_neighbour(constants->handle, setSize, reverse, (double *)mesh->x->data_d,
                   (double *)mesh->y->data_d, (double *)data->grx->data_d,
                   (double *)data->gsx->data_d, (double *)data->gry->data_d,
                   (double *)data->gsy->data_d);

  free(reverse);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(7, init_grad_args);
}
