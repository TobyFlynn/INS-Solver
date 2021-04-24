//
// auto-generated by op2.py
//

#include "cublas_v2.h"

#include  "op_lib_cpp.h"

//
// op_par_loop declarations
//
#ifdef OPENACC
#ifdef __cplusplus
extern "C" {
#endif
#endif
#ifdef OPENACC
#ifdef __cplusplus
}
#endif
#endif

#include "../blas_calls.h"

inline void cublas_poisson_bc2(cublasHandle_t handle, const int numCells,
                                const double *flux_d, double *rhs_d) {
  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 1.0;
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 21, &alpha, constants->gInterp_d, 15, flux_d, 21, &beta, rhs_d, 15);
}

void poisson_bc_blas2(INSData *data, Poisson_MF *poisson) {
  // Make sure OP2 data is in the right place
  op_arg poisson_args[] = {
    op_arg_dat(poisson->uFluxX, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->dudx, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 2, poisson_args);

  cublas_poisson_bc2(constants->handle, data->numCells, (double *)poisson->uFluxX->data_d,
                      (double *)poisson->dudx->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(2, poisson_args);
}
