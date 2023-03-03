#ifndef __INS_MPI_HELPER_FUNC_H
#define __INS_MPI_HELPER_FUNC_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

int compute_local_size(int global_size, int mpi_comm_size, int mpi_rank);

int compute_global_start(int global_size, int mpi_comm_size, int mpi_rank);

void scatter_DG_FP_array(DG_FP *g_array, DG_FP *l_array, int comm_size,
                          int g_size, int l_size, int elem_size);

void scatter_int_array(int *g_array, int *l_array, int comm_size, int g_size,
                       int l_size, int elem_size);

void gather_DG_FP_array(DG_FP *g_array, DG_FP *l_array, int comm_size,
                         int g_size, int l_size, int elem_size);

void gather_int_array(int *g_array, int *l_array, int comm_size, int g_size,
                      int l_size, int elem_size);

int get_global_mat_start_ind(int unknowns);

int get_global_element_start_ind(op_set set);

void gather_op2_DG_FP_array(DG_FP *g_array, DG_FP *l_array, int l_size,
                             int elem_size, int comm_size, int rank);

#endif
