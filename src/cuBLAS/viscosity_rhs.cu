#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_viscosity_rhs(cublasHandle_t handle, const int numCells, const double *mm_d,
                                 const double *qtt0_d, const double *qtt1_d,
                                 double *vis0_d, double *vis1_d) {
  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;

  for(int i = 0; i < numCells; i++) {
    const double *mm = mm_d + i * 15 * 15;
    const double *qtt0 = qtt0_d + i * 15;
    const double *qtt1 = qtt1_d + i * 15;
    double *vis0 = vis0_d + i * 15;
    double *vis1 = vis1_d + i * 15;

    cublasDgemv(handle, CUBLAS_OP_N, 15, 15, &alpha, mm, 15, qtt0, 1, &beta, vis0, 1);
    cublasDgemv(handle, CUBLAS_OP_N, 15, 15, &alpha, mm, 15, qtt1, 1, &beta, vis1, 1);
  }
}

void viscosity_rhs_blas(INSData *nsData, CubatureData *cubatureData) {
  // Make sure OP2 data is in the right place
  op_arg viscosity_rhs_args[] = {
    op_arg_dat(cubatureData->mm, -1, OP_ID, 15 * 15, "double", OP_READ),
    op_arg_dat(nsData->QTT[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->QTT[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->visRHS[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->visRHS[1], -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 5, viscosity_rhs_args);

  cublas_viscosity_rhs(constants->handle, nsData->numCells, (double *)cubatureData->mm->data_d, (double *)nsData->QTT[0]->data_d,
                      (double *)nsData->QTT[1]->data_d, (double *)nsData->visRHS[0]->data_d,
                      (double *)nsData->visRHS[1]->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(5, viscosity_rhs_args);
}
