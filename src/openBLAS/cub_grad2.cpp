#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_cub_grad2(const int numCells, const double *temp0, const double *temp1,
                              const double *temp2, const double *temp3, double *ux, double *uy) {
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 46, 1.0, constants->cubDr, 15, temp0, 46, 0.0, ux, 15);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 46, 1.0, constants->cubDs, 15, temp1, 46, 1.0, ux, 15);

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 46, 1.0, constants->cubDr, 15, temp2, 46, 0.0, uy, 15);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 46, 1.0, constants->cubDs, 15, temp3, 46, 1.0, uy, 15);
}

void cub_grad_blas2(INSData *data, CubatureData *cubatureData, op_dat ux, op_dat uy) {
  // Make sure OP2 data is in the right place
  op_arg grad_args[] = {
    op_arg_dat(cubatureData->op_temps[0], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[1], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[2], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[3], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(ux, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(uy, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 6, grad_args);

  openblas_cub_grad2(data->numCells, (double *)cubatureData->op_temps[0]->data,
                  (double *)cubatureData->op_temps[1]->data, (double *)cubatureData->op_temps[2]->data,
                  (double *)cubatureData->op_temps[3]->data, (double *)ux->data, (double *)uy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(6, grad_args);
}
