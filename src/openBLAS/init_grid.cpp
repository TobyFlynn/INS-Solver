#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_init_grid(const int numCells, const double *node_coords,
                               const int *cell2nodes, double *x, double *y,
                               double *xr, double *xs, double *yr, double *ys) {
  for(int c = 0; c < numCells; c++) {
    // Get nodes for this cell (on host)
    const double *n0 = &node_coords[2 * cell2nodes[3 * c]];
    const double *n1 = &node_coords[2 * cell2nodes[3 * c + 1]];
    const double *n2 = &node_coords[2 * cell2nodes[3 * c + 2]];

    double temp[15];
    double *x_c = x + c * 15;
    double *y_c = y + c * 15;
    double *xr_c = xr + c * 15;
    double *xs_c = xs + c * 15;
    double *yr_c = yr + c * 15;
    double *ys_c = ys + c * 15;

    cblas_dcopy(15, ones, 1, x_c, 1);
    cblas_daxpy(15, 1.0, r, 1, x_c, 1);
    cblas_dscal(15, 0.5 * n1[0], x_c, 1);
    cblas_dcopy(15, ones, 1, temp, 1);
    cblas_daxpy(15, 1.0, s, 1, temp, 1);
    cblas_daxpy(15, 0.5 * n2[0], temp, 1, x_c, 1);
    cblas_dcopy(15, s, 1, temp, 1);
    cblas_daxpy(15, 1.0, r, 1, temp, 1);
    cblas_daxpy(15, -0.5 * n0[0], temp, 1, x_c, 1);

    cblas_dcopy(15, ones, 1, y_c, 1);
    cblas_daxpy(15, 1.0, r, 1, y_c, 1);
    cblas_dscal(15, 0.5 * n1[1], y_c, 1);
    cblas_dcopy(15, ones, 1, temp, 1);
    cblas_daxpy(15, 1.0, s, 1, temp, 1);
    cblas_daxpy(15, 0.5 * n2[1], temp, 1, y_c, 1);
    cblas_dcopy(15, s, 1, temp, 1);
    cblas_daxpy(15, 1.0, r, 1, temp, 1);
    cblas_daxpy(15, -0.5 * n0[1], temp, 1, y_c, 1);

    // xr = Dr * x
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dr, 15, x_c, 1, 0.0, xr_c, 1);
    // xs = Ds * x
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Ds, 15, x_c, 1, 0.0, xs_c, 1);
    // yr = Dr * y
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dr, 15, y_c, 1, 0.0, yr_c, 1);
    // ys = Ds * y
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Ds, 15, y_c, 1, 0.0, ys_c, 1);
  }
}

void init_grid_blas(INSData *nsData) {
  // Make sure OP2 data is in the right place
  op_arg init_grid_args[] = {
    op_arg_dat(nsData->x, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->y, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->xr, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->xs, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->yr, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->ys, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 6, init_grid_args);

  openblas_init_grid(nsData->numCells, (double *)nsData->node_coords->data,
                     (int *)nsData->cell2nodes->map, (double *)nsData->x->data,
                     (double *)nsData->y->data, (double *)nsData->xr->data,
                     (double *)nsData->xs->data, (double *)nsData->yr->data,
                     (double *)nsData->ys->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(6, init_grid_args);
}
