#include "mpi.h"

#include <memory>

int compute_local_size(int global_size, int mpi_comm_size, int mpi_rank) {
  int local_size = global_size / mpi_comm_size;
  int remainder = global_size % mpi_comm_size;

  if (mpi_rank < remainder) {
    local_size = local_size + 1;
  }
  return local_size;
}

void scatter_double_array(double *g_array, double *l_array, int comm_size,
                          int g_size, int l_size, int elem_size) {
  int *sendcnts = (int *)malloc(comm_size * sizeof(int));
  int *displs = (int *)malloc(comm_size * sizeof(int));
  int disp = 0;

  for (int i = 0; i < comm_size; i++) {
    sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
  }
  for (int i = 0; i < comm_size; i++) {
    displs[i] = disp;
    disp = disp + sendcnts[i];
  }

  MPI_Scatterv(g_array, sendcnts, displs, MPI_DOUBLE, l_array,
               l_size * elem_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(sendcnts);
  free(displs);
}

void scatter_int_array(int *g_array, int *l_array, int comm_size, int g_size,
                       int l_size, int elem_size) {
  int *sendcnts = (int *)malloc(comm_size * sizeof(int));
  int *displs = (int *)malloc(comm_size * sizeof(int));
  int disp = 0;

  for (int i = 0; i < comm_size; i++) {
    sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
  }
  for (int i = 0; i < comm_size; i++) {
    displs[i] = disp;
    disp = disp + sendcnts[i];
  }

  MPI_Scatterv(g_array, sendcnts, displs, MPI_INT, l_array, l_size * elem_size,
               MPI_INT, 0, MPI_COMM_WORLD);

  free(sendcnts);
  free(displs);
}

void gather_double_array(double *g_array, double *l_array, int comm_size,
                         int g_size, int l_size, int elem_size) {
  int *sendcnts = (int *)malloc(comm_size * sizeof(int));
  int *displs = (int *)malloc(comm_size * sizeof(int));
  int disp = 0;

  for (int i = 0; i < comm_size; i++) {
    sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
  }
  for (int i = 0; i < comm_size; i++) {
    displs[i] = disp;
    disp = disp + sendcnts[i];
  }

  MPI_Gatherv(l_array, l_size * elem_size, MPI_DOUBLE, g_array, sendcnts,
              displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(sendcnts);
  free(displs);
}
