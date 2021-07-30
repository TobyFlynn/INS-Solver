#include "cblas.h"

#include "op_seq.h"
#include "blas_calls.h"

inline void openblas_init_gauss_grad_neighbour(const int numCells, const int *reverse, const double *x,
                                const double *y, double *gxr, double *gxs,
                                double *gyr, double *gys) {
  for(int c = 0; c < numCells; c++) {
    const double *x_c = x + c * 6;
    const double *y_c = y + c * 6;
    double *gxr_c = gxr + c * 12;
    double *gxs_c = gxs + c * 12;
    double *gyr_c = gyr + c * 12;
    double *gys_c = gys + c * 12;

    if(reverse[c * 3]) {
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF0DrR, 6, x_c, 1, 0.0, gxr_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF0DsR, 6, x_c, 1, 0.0, gxs_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF0DrR, 6, y_c, 1, 0.0, gyr_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF0DsR, 6, y_c, 1, 0.0, gys_c, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF0Dr, 6, x_c, 1, 0.0, gxr_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF0Ds, 6, x_c, 1, 0.0, gxs_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF0Dr, 6, y_c, 1, 0.0, gyr_c, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF0Ds, 6, y_c, 1, 0.0, gys_c, 1);
    }

    if(reverse[c * 3 + 1]) {
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF1DrR, 6, x_c, 1, 0.0, gxr_c + 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF1DsR, 6, x_c, 1, 0.0, gxs_c + 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF1DrR, 6, y_c, 1, 0.0, gyr_c + 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF1DsR, 6, y_c, 1, 0.0, gys_c + 4, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF1Dr, 6, x_c, 1, 0.0, gxr_c + 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF1Ds, 6, x_c, 1, 0.0, gxs_c + 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF1Dr, 6, y_c, 1, 0.0, gyr_c + 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF1Ds, 6, y_c, 1, 0.0, gys_c + 4, 1);
    }

    if(reverse[c * 3 + 2]) {
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF2DrR, 6, x_c, 1, 0.0, gxr_c + 2 * 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF2DsR, 6, x_c, 1, 0.0, gxs_c + 2 * 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF2DrR, 6, y_c, 1, 0.0, gyr_c + 2 * 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF2DsR, 6, y_c, 1, 0.0, gys_c + 2 * 4, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF2Dr, 6, x_c, 1, 0.0, gxr_c + 2 * 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF2Ds, 6, x_c, 1, 0.0, gxs_c + 2 * 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF2Dr, 6, y_c, 1, 0.0, gyr_c + 2 * 4, 1);
      cblas_dgemv(CblasColMajor, CblasTrans, 6, 4, 1.0, constants->gF2Ds, 6, y_c, 1, 0.0, gys_c + 2 * 4, 1);
    }
  }
}

void init_gauss_grad_neighbour_blas(DGMesh *mesh, INSData *data) {
  // Make sure OP2 data is in the right place
  op_arg init_grad_args[] = {
    op_arg_dat(mesh->x, -1, OP_ID, 6, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, 6, "double", OP_READ),
    op_arg_dat(data->reverse, -1, OP_ID, 3, "int", OP_READ),
    op_arg_dat(data->grx, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(data->gsx, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(data->gry, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(data->gsy, -1, OP_ID, 21, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 7, init_grad_args);

  int setSize = mesh->x->set->size;

  openblas_init_gauss_grad_neighbour(setSize, (int *)data->reverse->data, (double *)mesh->x->data,
                   (double *)mesh->y->data, (double *)data->grx->data,
                   (double *)data->gsx->data, (double *)data->gry->data,
                   (double *)data->gsy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(7, init_grad_args);
}
