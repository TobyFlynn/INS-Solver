#include "matrices/2d/poisson_matrix_2d.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

#include "dg_utils.h"
#include "timing.h"

extern Timing *timer;

void PoissonMatrix2D::set_glb_ind() {
  timer->startTimer("PoissonMat - glb_ind");
  int global_ind = 0;
  #ifdef INS_MPI
  timer->startTimer("PoissonMat - unknowns");
  op_arg op2_args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->order->set, 1, op2_args);
  const int setSize_ = mesh->order->set->size;
  int *tempOrder_ = (int *)malloc(setSize_ * sizeof(int));
  cudaMemcpy(tempOrder_, mesh->order->data_d, setSize_ * sizeof(int), cudaMemcpyDeviceToHost);
  unknowns = 0;
  for(int i = 0; i < setSize_; i++) {
    int Np, Nfp;
    DGUtils::numNodes2D(tempOrder_[i], &Np, &Nfp);
    unknowns += Np;
  }
  free(tempOrder_);
  op_mpi_set_dirtybit_cuda(1, op2_args);
  global_ind = get_global_mat_start_ind(unknowns);
  #endif
  op_arg args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 2, args);

  const int setSize = mesh->cells->size;
  int *tempOrder = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(tempOrder, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  int *data_ptr = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(data_ptr, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  int ind = global_ind;
  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::numNodes2D(tempOrder[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  cudaMemcpy(glb_ind->data_d, data_ptr, setSize * sizeof(int), cudaMemcpyHostToDevice);

  op_mpi_set_dirtybit_cuda(2, args);
  free(data_ptr);
  free(tempOrder);
  timer->endTimer("PoissonMat - glb_ind");
}
