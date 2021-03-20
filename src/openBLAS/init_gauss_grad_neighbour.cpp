#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_init_gauss_grad_neighbour(const int numCells, const int *reverse, const double *x,
                                const double *y, double *gxr, double *gxs,
                                double *gyr, double *gys) {
  for(int c = 0; c < numCells; c++) {
    const double *x_c = x + c * 15;
    const double *y_c = y + c * 15;
    double *gxr_c = gxr + c * 21;
    double *gxs_c = gxs + c * 21;
    double *gyr_c = gyr + c * 21;
    double *gys_c = gys + c * 21;

    if(reverse[c * 3]) {
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF0DrR, 15, x_c, 1, 0.0, gxr_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF0DsR, 15, x_c, 1, 0.0, gxs_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF0DrR, 15, y_c, 1, 0.0, gyr_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF0DsR, 15, y_c, 1, 0.0, gys_c, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF0Dr, 15, x_c, 1, 0.0, gxr_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF0Ds, 15, x_c, 1, 0.0, gxs_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF0Dr, 15, y_c, 1, 0.0, gyr_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF0Ds, 15, y_c, 1, 0.0, gys_c, 1);
    }

    if(reverse[c * 3 + 1]) {
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF1DrR, 15, x_c, 1, 0.0, gxr_c + 7, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF1DsR, 15, x_c, 1, 0.0, gxs_c + 7, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF1DrR, 15, y_c, 1, 0.0, gyr_c + 7, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF1DsR, 15, y_c, 1, 0.0, gys_c + 7, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF1Dr, 15, x_c, 1, 0.0, gxr_c + 7, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF1Ds, 15, x_c, 1, 0.0, gxs_c + 7, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF1Dr, 15, y_c, 1, 0.0, gyr_c + 7, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF1Ds, 15, y_c, 1, 0.0, gys_c + 7, 1);
    }

    if(reverse[c * 3 + 2]) {
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF2DrR, 15, x_c, 1, 0.0, gxr_c + 14, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF2DsR, 15, x_c, 1, 0.0, gxs_c + 14, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF2DrR, 15, y_c, 1, 0.0, gyr_c + 14, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF2DsR, 15, y_c, 1, 0.0, gys_c + 14, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF2Dr, 15, x_c, 1, 0.0, gxr_c + 14, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF2Ds, 15, x_c, 1, 0.0, gxs_c + 14, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF2Dr, 15, y_c, 1, 0.0, gyr_c + 14, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 15, 7, 1.0, gF2Ds, 15, y_c, 1, 0.0, gys_c + 14, 1);
    }
  }
}

void init_gauss_grad_neighbour_blas(INSData *nsData, GaussData *gaussData) {
  // Make sure OP2 data is in the right place
  op_arg init_grad_args[] = {
    op_arg_dat(nsData->x, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->y, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(gaussData->rx, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->sx, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->ry, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->sy, -1, OP_ID, 21, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 6, init_grad_args);

  int *reverse = (int *)malloc(3 * op_get_size(nsData->cells) * sizeof(int));
  op_fetch_data(gaussData->reverse, reverse);

  openblas_init_gauss_grad_neighbour(nsData->numCells, reverse, (double *)nsData->x->data,
                   (double *)nsData->y->data, (double *)gaussData->rx->data,
                   (double *)gaussData->sx->data, (double *)gaussData->ry->data,
                   (double *)gaussData->sy->data);

  free(reverse);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(6, init_grad_args);
}