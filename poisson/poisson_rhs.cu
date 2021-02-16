#include "poisson_rhs.h"

#include "op_seq.h"

#include "operators.h"
#include "blas_calls.h"

#include "kernels/pRHS_faces.h"
#include "kernels/pRHS_bc.h"
#include "kernels/pRHS_du.h"
#include "kernels/pRHS_fluxq.h"
#include "kernels/pRHS_J.h"
#include "kernels/pRHS_qbc.h"

void poisson_rhs(const double *u, double *rhs) {
  op_arg u_copy_args[] = {
    op_arg_dat(data->pU, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, u_copy_args);
  cudaMemcpy(data->pU->data_d, u, data->numCells * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, u_copy_args);

  op_par_loop(pRHS_faces, "pRHS_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->pU, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->pExRHS[0], -2, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(pRHS_bc, "pRHS_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->pU, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->pExRHS[0], 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(pRHS_du, "pRHS_du", data->cells,
              op_arg_dat(data->nx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ny, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->fscale, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pU, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pExRHS[0], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->pDu, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->pFluxXu, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(data->pFluxYu, -1, OP_ID, 15, "double", OP_WRITE));

  grad(data, data->pU, data->pDuDx, data->pDuDy);

  // qx and qy stored in pDuDx and pDuDy
  poisson_rhs_blas1(data);

  op_par_loop(pRHS_faces, "pRHS_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->pDuDx, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->pExRHS[0], -2, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(pRHS_faces, "pRHS_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->pDuDy, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->pExRHS[1], -2, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(pRHS_qbc, "pRHS_qbc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->pDuDx, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->pExRHS[0], 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(pRHS_qbc, "pRHS_qbc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->pDuDy, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->pExRHS[1], 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(pRHS_fluxq, "pRHS_fluxq", data->cells,
              op_arg_dat(data->nx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ny, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->fscale, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pTau, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pDu, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pDuDx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pDuDy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pExRHS[0], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->pExRHS[1], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(data->pFluxQ, -1, OP_ID, 15, "double", OP_WRITE));

  div(data, data->pDuDx, data->pDuDy, data->pDivQ);

  poisson_rhs_blas2(data);

  op_par_loop(pRHS_J, "pRHS_J", data->cells,
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pRHSU, -1, OP_ID, 15, "double", OP_RW));

  op_arg rhs_copy_args[] = {
    op_arg_dat(data->pRHSU, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, rhs_copy_args);
  cudaMemcpy(rhs, data->pRHSU->data_d, data->numCells * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, rhs_copy_args);
}
