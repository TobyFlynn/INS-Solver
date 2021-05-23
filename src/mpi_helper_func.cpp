#include "op_seq.h"

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

int compute_global_start(int global_size, int mpi_comm_size, int mpi_rank) {
  int start = 0;
  for (int i = 0; i < mpi_rank; i++) {
    start += compute_local_size(global_size, mpi_comm_size, i);
  }
  return start;
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

int get_global_start_index(op_set set) {
  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int *sizes = (int *)malloc(comm_size * sizeof(int));

  MPI_Allgather(&set->size, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);

  int index = 0;
  for(int i = 0; i < rank; i++) {
    index += sizes[i];
  }

  free(sizes);
  return index;
}

void gather_op2_double_array(double *g_array, double *l_array, int l_size,
                             int elem_size, int comm_size, int rank) {
  int *sizes = (int *)malloc(comm_size * sizeof(int));
  MPI_Allgather(&l_size, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);

  int *sendcnts = (int *)malloc(comm_size * sizeof(int));
  int *displs = (int *)malloc(comm_size * sizeof(int));
  int disp = 0;

  for (int i = 0; i < comm_size; i++) {
    sendcnts[i] = elem_size * sizes[i];
  }
  for (int i = 0; i < comm_size; i++) {
    displs[i] = disp;
    disp = disp + sendcnts[i];
  }

  MPI_Gatherv(l_array, sizes[rank] * elem_size, MPI_DOUBLE, g_array, sendcnts,
              displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(sendcnts);
  free(displs);
  free(sizes);
}
