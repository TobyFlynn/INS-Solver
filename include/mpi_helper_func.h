#ifndef __INS_MPI_HELPER_FUNC_H
#define __INS_MPI_HELPER_FUNC_H

#include "op_seq.h"

int compute_local_size(int global_size, int mpi_comm_size, int mpi_rank);

int compute_global_start(int global_size, int mpi_comm_size, int mpi_rank);

void scatter_double_array(double *g_array, double *l_array, int comm_size,
                          int g_size, int l_size, int elem_size);

void scatter_int_array(int *g_array, int *l_array, int comm_size, int g_size,
                       int l_size, int elem_size);

void gather_double_array(double *g_array, double *l_array, int comm_size,
                         int g_size, int l_size, int elem_size);

int get_global_start_index(op_set set);

#endif
