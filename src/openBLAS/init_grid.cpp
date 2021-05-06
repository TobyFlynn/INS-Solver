#include "cblas.h"

#include "op_seq.h"
#include "blas_calls.h"

inline void openblas_init_grid(const int numCells, const double *nodeX,
                               const double *nodeY, double *x, double *y,
                               double *xr, double *xs, double *yr, double *ys) {
  for(int c = 0; c < numCells; c++) {
    // Get nodes for this cell (on host)
    const double n0[] = {nodeX[c * 3], nodeY[3 * c]};
    const double n1[] = {nodeX[c * 3 + 1], nodeY[3 * c + 1]};
    const double n2[] = {nodeX[c * 3 + 2], nodeY[3 * c + 2]};

    double temp[15];
    double *x_c = x + c * 15;
    double *y_c = y + c * 15;
    double *xr_c = xr + c * 15;
    double *xs_c = xs + c * 15;
    double *yr_c = yr + c * 15;
    double *ys_c = ys + c * 15;

    cblas_dcopy(15, constants->ones, 1, x_c, 1);
    cblas_daxpy(15, 1.0, constants->r, 1, x_c, 1);
    cblas_dscal(15, 0.5 * n1[0], x_c, 1);
    cblas_dcopy(15, constants->ones, 1, temp, 1);
    cblas_daxpy(15, 1.0, constants->s, 1, temp, 1);
    cblas_daxpy(15, 0.5 * n2[0], temp, 1, x_c, 1);
    cblas_dcopy(15, constants->s, 1, temp, 1);
    cblas_daxpy(15, 1.0, constants->r, 1, temp, 1);
    cblas_daxpy(15, -0.5 * n0[0], temp, 1, x_c, 1);

    cblas_dcopy(15, constants->ones, 1, y_c, 1);
    cblas_daxpy(15, 1.0, constants->r, 1, y_c, 1);
    cblas_dscal(15, 0.5 * n1[1], y_c, 1);
    cblas_dcopy(15, constants->ones, 1, temp, 1);
    cblas_daxpy(15, 1.0, constants->s, 1, temp, 1);
    cblas_daxpy(15, 0.5 * n2[1], temp, 1, y_c, 1);
    cblas_dcopy(15, constants->s, 1, temp, 1);
    cblas_daxpy(15, 1.0, constants->r, 1, temp, 1);
    cblas_daxpy(15, -0.5 * n0[1], temp, 1, y_c, 1);
  }

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->Dr, 15, x, 15, 0.0, xr, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->Ds, 15, x, 15, 0.0, xs, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->Dr, 15, y, 15, 0.0, yr, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->Ds, 15, y, 15, 0.0, ys, 15);
}

void init_grid_blas(INSData *nsData) {
  // Make sure OP2 data is in the right place
  op_arg init_grid_args[] = {
    op_arg_dat(nsData->nodeX, -1, OP_ID, 3, "double", OP_READ),
    op_arg_dat(nsData->nodeY, -1, OP_ID, 3, "double", OP_READ),
    op_arg_dat(nsData->x, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->y, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->rx, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->sx, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->ry, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->sy, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 8, init_grid_args);

  int setSize = nsData->x->set->size;

  openblas_init_grid(setSize, (double *)nsData->nodeX->data,
                     (double *)nsData->nodeY->data, (double *)nsData->x->data,
                     (double *)nsData->y->data, (double *)nsData->rx->data,
                     (double *)nsData->sx->data, (double *)nsData->ry->data,
                     (double *)nsData->sy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(8, init_grid_args);
}
