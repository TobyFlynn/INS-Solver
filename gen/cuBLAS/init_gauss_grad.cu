#include "cublas_v2.h"

#include "op_seq.h"
#include "blas_calls.h"

inline void cublas_init_gauss_grad(cublasHandle_t handle, const int numCells,
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
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0Dr_d, 6, x, 1, &beta, gxr, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0Ds_d, 6, x, 1, &beta, gxs, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0Dr_d, 6, y, 1, &beta, gyr, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF0Ds_d, 6, y, 1, &beta, gys, 1);

    // Face 1
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1Dr_d, 6, x, 1, &beta, gxr + 6, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1Ds_d, 6, x, 1, &beta, gxs + 6, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1Dr_d, 6, y, 1, &beta, gyr + 6, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF1Ds_d, 6, y, 1, &beta, gys + 6, 1);

    // Face 2
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2Dr_d, 6, x, 1, &beta, gxr + 2 * 6, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2Ds_d, 6, x, 1, &beta, gxs + 2 * 6, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2Dr_d, 6, y, 1, &beta, gyr + 2 * 6, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 6, 10, &alpha, constants->gF2Ds_d, 6, y, 1, &beta, gys + 2 * 6, 1);
  }
}

void init_gauss_grad_blas(DGMesh *mesh, INSData *data) {
  // Make sure OP2 data is in the right place
  op_arg init_grad_args[] = {
    op_arg_dat(mesh->x, -1, OP_ID, 10, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, 10, "double", OP_READ),
    op_arg_dat(data->grx, -1, OP_ID, 18, "double", OP_WRITE),
    op_arg_dat(data->gsx, -1, OP_ID, 18, "double", OP_WRITE),
    op_arg_dat(data->gry, -1, OP_ID, 18, "double", OP_WRITE),
    op_arg_dat(data->gsy, -1, OP_ID, 18, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 6, init_grad_args);

  int setSize = mesh->x->set->size;

  cublas_init_gauss_grad(constants->handle, setSize, (double *)mesh->x->data_d,
                   (double *)mesh->y->data_d, (double *)data->grx->data_d,
                   (double *)data->gsx->data_d, (double *)data->gry->data_d,
                   (double *)data->gsy->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(6, init_grad_args);
}
