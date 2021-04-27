#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

// inline void cublas_poisson_rhs1(cublasHandle_t handle, const int numCells,
//                                 const double *fluxX_d, const double *fluxY_d,
//                                 double *dudx_d, double *dudy_d, double *qx_d,
//                                 double *qy_d, double *gradx_d, double *grady_d) {
//   // CUBLAS_OP_T because cublas is column major but constants are stored row major
//   double alpha = -1.0;
//   double alpha2 = 1.0;
//   double beta = 1.0;
//   double beta2 = 0.0;
//   cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, constants->invMass_d, 15, dudx_d, 15, &beta2, gradx_d, 15);
//   cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, constants->invMass_d, 15, dudy_d, 15, &beta2, grady_d, 15);
//
//   cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 21, &alpha, constants->gInterp_d, 15, fluxX_d, 21, &beta, dudx_d, 15);
//   cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 21, &alpha, constants->gInterp_d, 15, fluxY_d, 21, &beta, dudy_d, 15);
//
//   cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, constants->invMass_d, 15, dudx_d, 15, &beta2, qx_d, 15);
//   cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha2, constants->invMass_d, 15, dudy_d, 15, &beta2, qy_d, 15);
// }

void poisson_mf2_blas(INSData *data, Poisson_MF2 *poisson) {
  // Make sure OP2 data is in the right place
  op_arg poisson_args[] = {
    op_arg_dat(poisson->u, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(poisson->op1, -1, OP_ID, 15 * 15, "double", OP_READ),
    op_arg_dat(poisson->rhs, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 3, poisson_args);

  // cublas_poisson_rhs1(constants->handle, data->numCells, (double *)poisson->u->data_d,
  //                     (double *)poisson->op1->data_d, (double *)poisson->rhs->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(3, poisson_args);
}
