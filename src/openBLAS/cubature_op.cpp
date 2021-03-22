#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_cubature_op(const int numCells, const double *Dx, const double *Dy,
                                 const double *temp, const double *temp2, double *OP) {
  for(int i = 0; i < numCells; i++) {
    const double *Dx_c = Dx + i * 46 * 15;
    const double *Dy_c = Dy + i * 46 * 15;
    const double *temp_c = temp + i * 46 * 15;
    const double *temp2_c = temp2 + i * 46 * 15;
    double *OP_c = OP + i * 15 * 15;

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 15, 15, 46, 1.0, Dx_c, 15, temp_c, 15, 0.0, OP_c, 15);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 15, 15, 46, 1.0, Dy_c, 15, temp2_c, 15, 1.0, OP_c, 15);
  }
}

void cubature_op_blas(INSData *nsData, CubatureData *cubData) {
  // Make sure OP2 data is in the right place
  op_arg op_cubature_args[] = {
    op_arg_dat(cubData->Dx, -1, OP_ID, 46 * 15, "double", OP_READ),
    op_arg_dat(cubData->Dy, -1, OP_ID, 46 * 15, "double", OP_READ),
    op_arg_dat(cubData->temp, -1, OP_ID, 46 * 15, "double", OP_READ),
    op_arg_dat(cubData->temp2, -1, OP_ID, 46 * 15, "double", OP_READ),
    op_arg_dat(cubData->OP, -1, OP_ID, 15 * 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 5, op_cubature_args);

  openblas_cubature_op(nsData->numCells, (double *)cubData->Dx->data,
                     (double *)cubData->Dy->data, (double *)cubData->temp->data,
                     (double *)cubData->temp2->data, (double *)cubData->OP->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(5, op_cubature_args);
}
