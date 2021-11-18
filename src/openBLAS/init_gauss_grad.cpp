#include "cblas.h"

#include "op_seq.h"
#include "blas_calls.h"

inline void openblas_init_gauss_grad(const int numCells, const double *x,
                                const double *y, double *gxr, double *gxs,
                                double *gyr, double *gys) {
  for(int c = 0; c < numCells; c++) {
    const double *x_c = x + c * 15;
    const double *y_c = y + c * 15;
    double *gxr_c = gxr + c * 21;
    double *gxs_c = gxs + c * 21;
    double *gyr_c = gyr + c * 21;
    double *gys_c = gys + c * 21;

    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF0Dr, 7, x_c, 1, 0.0, gxr_c, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF0Ds, 7, x_c, 1, 0.0, gxs_c, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF0Dr, 7, y_c, 1, 0.0, gyr_c, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF0Ds, 7, y_c, 1, 0.0, gys_c, 1);

    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF1Dr, 7, x_c, 1, 0.0, gxr_c + 7, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF1Ds, 7, x_c, 1, 0.0, gxs_c + 7, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF1Dr, 7, y_c, 1, 0.0, gyr_c + 7, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF1Ds, 7, y_c, 1, 0.0, gys_c + 7, 1);

    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF2Dr, 7, x_c, 1, 0.0, gxr_c + 14, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF2Ds, 7, x_c, 1, 0.0, gxs_c + 14, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF2Dr, 7, y_c, 1, 0.0, gyr_c + 14, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 7, 15, 1.0, constants->gF2Ds, 7, y_c, 1, 0.0, gys_c + 14, 1);
  }
}

void init_gauss_grad_blas(DGMesh *mesh, GaussData *gaussData) {
  // Make sure OP2 data is in the right place
  op_arg init_grad_args[] = {
    op_arg_dat(mesh->x, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(gaussData->rx, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->sx, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->ry, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->sy, -1, OP_ID, 21, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 6, init_grad_args);

  int setSize = mesh->x->set->size;

  openblas_init_gauss_grad(setSize, (double *)mesh->x->data,
                   (double *)mesh->y->data, (double *)gaussData->rx->data,
                   (double *)gaussData->sx->data, (double *)gaussData->ry->data,
                   (double *)gaussData->sy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(6, init_grad_args);
}
