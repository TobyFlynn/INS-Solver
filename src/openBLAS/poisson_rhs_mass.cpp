#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_rhs_mass(const int numCells, const double *mm,
                                      const double *u, double *rhs, double factor) {
  for(int i = 0; i < numCells; i++) {
    const double *mm_c = mm + i * 15 * 15;
    const double *u_c = u + i * 15;
    double *rhs_c = rhs + i * 15;

    cblas_dgemv(CblasColMajor, CblasNoTrans, 15, 15, factor, mm_c, 15, u_c, 1, 1.0, rhs_c, 1);
  }
}

void poisson_rhs_mass_blas(INSData *data, CubatureData *cubatureData, Poisson_MF *poisson, double factor) {
  // Make sure OP2 data is in the right place
  op_arg mass_args[] = {
    op_arg_dat(cubatureData->mm, -1, OP_ID, 15 * 15, "double", OP_READ),
    op_arg_dat(poisson->u, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(poisson->rhs, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges(data->cells, 3, mass_args);

  openblas_poisson_rhs_mass(data->numCells, (double *)cubatureData->mm->data,
                          (double *)poisson->u->data, (double *)poisson->rhs->data, factor);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(3, mass_args);
}
