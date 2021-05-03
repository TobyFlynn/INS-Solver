#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_viscosity_rhs(const int numCells, const double *mm, const double *qtt0,
                                   const double *qtt1, double *vis0, double *vis1) {
  for(int i = 0; i < numCells; i++) {
    const double *mm_c = mm + i * 15 * 15;
    const double *qtt0_c = qtt0 + i * 15;
    const double *qtt1_c = qtt1 + i * 15;
    double *vis0_c = vis0 + i * 15;
    double *vis1_c = vis1 + i * 15;

    cblas_dgemv(CblasColMajor, CblasNoTrans, 15, 15, 1.0, mm_c, 15, qtt0_c, 1, 0.0, vis0_c, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 15, 15, 1.0, mm_c, 15, qtt1_c, 1, 0.0, vis1_c, 1);
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
  op_mpi_halo_exchanges(nsData->cells, 5, viscosity_rhs_args);

  openblas_viscosity_rhs(nsData->numCells, (double *)cubatureData->mm->data, (double *)nsData->QTT[0]->data,
                      (double *)nsData->QTT[1]->data, (double *)nsData->visRHS[0]->data,
                      (double *)nsData->visRHS[1]->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(5, viscosity_rhs_args);
}
