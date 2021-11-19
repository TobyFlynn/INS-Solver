#include "cblas.h"

#include "op_seq.h"
#include "blas_calls.h"

inline void openblas_init_gauss_grad_neighbour(const int numCells, const int *reverse, const double *x,
                                const double *y, double *gxr, double *gxs,
                                double *gyr, double *gys) {
  for(int c = 0; c < numCells; c++) {
    const double *x_c = x + c * DG_NP;
    const double *y_c = y + c * DG_NP;
    double *gxr_c = gxr + c * DG_G_NP;
    double *gxs_c = gxs + c * DG_G_NP;
    double *gyr_c = gyr + c * DG_G_NP;
    double *gys_c = gys + c * DG_G_NP;

    if(reverse[c * 3]) {
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF0DrR, DG_GF_NP, x_c, 1, 0.0, gxr_c, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF0DsR, DG_GF_NP, x_c, 1, 0.0, gxs_c, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF0DrR, DG_GF_NP, y_c, 1, 0.0, gyr_c, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF0DsR, DG_GF_NP, y_c, 1, 0.0, gys_c, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF0Dr, DG_GF_NP, x_c, 1, 0.0, gxr_c, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF0Ds, DG_GF_NP, x_c, 1, 0.0, gxs_c, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF0Dr, DG_GF_NP, y_c, 1, 0.0, gyr_c, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF0Ds, DG_GF_NP, y_c, 1, 0.0, gys_c, 1);
    }

    if(reverse[c * 3 + 1]) {
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF1DrR, DG_GF_NP, x_c, 1, 0.0, gxr_c + DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF1DsR, DG_GF_NP, x_c, 1, 0.0, gxs_c + DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF1DrR, DG_GF_NP, y_c, 1, 0.0, gyr_c + DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF1DsR, DG_GF_NP, y_c, 1, 0.0, gys_c + DG_GF_NP, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF1Dr, DG_GF_NP, x_c, 1, 0.0, gxr_c + DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF1Ds, DG_GF_NP, x_c, 1, 0.0, gxs_c + DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF1Dr, DG_GF_NP, y_c, 1, 0.0, gyr_c + DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF1Ds, DG_GF_NP, y_c, 1, 0.0, gys_c + DG_GF_NP, 1);
    }

    if(reverse[c * 3 + 2]) {
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF2DrR, DG_GF_NP, x_c, 1, 0.0, gxr_c + 2 * DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF2DsR, DG_GF_NP, x_c, 1, 0.0, gxs_c + 2 * DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF2DrR, DG_GF_NP, y_c, 1, 0.0, gyr_c + 2 * DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF2DsR, DG_GF_NP, y_c, 1, 0.0, gys_c + 2 * DG_GF_NP, 1);
    } else {
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF2Dr, DG_GF_NP, x_c, 1, 0.0, gxr_c + 2 * DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF2Ds, DG_GF_NP, x_c, 1, 0.0, gxs_c + 2 * DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF2Dr, DG_GF_NP, y_c, 1, 0.0, gyr_c + 2 * DG_GF_NP, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans, DG_GF_NP, DG_NP, 1.0, constants->gF2Ds, DG_GF_NP, y_c, 1, 0.0, gys_c + 2 * DG_GF_NP, 1);
    }
  }
}

void init_gauss_grad_neighbour_blas(DGMesh *mesh, GaussData *gaussData) {
  // Make sure OP2 data is in the right place
  op_arg init_grad_args[] = {
    op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(gaussData->reverse, -1, OP_ID, 3, "int", OP_READ),
    op_arg_dat(gaussData->rx, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
    op_arg_dat(gaussData->sx, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
    op_arg_dat(gaussData->ry, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
    op_arg_dat(gaussData->sy, -1, OP_ID, DG_G_NP, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 7, init_grad_args);

  int setSize = mesh->x->set->size;

  openblas_init_gauss_grad_neighbour(setSize, (int *)gaussData->reverse->data, (double *)mesh->x->data,
                   (double *)mesh->y->data, (double *)gaussData->rx->data,
                   (double *)gaussData->sx->data, (double *)gaussData->ry->data,
                   (double *)gaussData->sy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(7, init_grad_args);
}
