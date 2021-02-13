#include "poisson_rhs.h"

#include "op_seq.h"

#include "kernels/pRHS_faces.h"
#include "kernels/pRHS_bc.h"
#include "kernels/pRHS_du.h"

void poisson_rhs(const double *u, double *rhs) {
  op_arg u_copy_args[] = {
    op_arg_dat(data->pU, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 1, u_copy_args);
  memcpy(data->pU->data, u, data->numCells * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, u_copy_args);

  op_par_loop(pRHS_faces, "pRHS_faces", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->pU, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(data->pExU, -2, data->edge2cells, 15, "double", OP_INC));

  op_par_loop(pRHS_bc, "pRHS_bc", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->pU, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->pExU, 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(pRHS_du, "pRHS_du", data->cells,
              op_arg_dat(data->pU, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pExU, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->pDu, -1, OP_ID, 15, "double", OP_WRITE);
}
