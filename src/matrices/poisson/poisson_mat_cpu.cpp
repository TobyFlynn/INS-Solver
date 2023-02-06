#include "poisson_matrix.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

#include "dg_utils.h"
#include "timing.h"

extern Timing *timer;

int PoissonMatrix2D::get_local_unknowns() {
  timer->startTimer("PoissonMat - unknowns");
  op_arg op2_args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(mesh->order->set, 1, op2_args);
  const int setSize = mesh->order->set->size;
  const int *p = (int *)mesh->order->data;
  int local_unkowns = 0;
  for(int i = 0; i < setSize; i++) {
    int Np, Nfp;
    DGUtils::numNodes2D(p[i], &Np, &Nfp);
    local_unkowns += Np;
  }
  op_mpi_set_dirtybit(1, op2_args);
  timer->endTimer("PoissonMat - unknowns");
  return local_unkowns;
}

void PoissonMatrix2D::set_glb_ind() {
  timer->startTimer("PoissonMat - glb_ind");
  int global_ind = 0;
  #ifdef INS_MPI
  global_ind = get_global_mat_start_ind(unknowns);
  #endif
  op_arg args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges(mesh->cells, 2, args);

  const int *p = (int *)mesh->order->data;
  int *data_ptr = (int *)glb_ind->data;
  int ind = global_ind;
  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::numNodes2D(p[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  op_mpi_set_dirtybit(2, args);
  timer->endTimer("PoissonMat - glb_ind");
}
