#include "hypre_utils.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"

int get_global_start_ind(const int local_unknowns) {
  return get_global_mat_start_ind(local_unknowns);
}
#else
int get_global_start_ind(const int local_unknowns) {
  return 0;
}
#endif

// Get an OP2 dat as a HYPRE vector
void HYPREUtils::dat_to_new_vec(op_dat v_dat, HYPRE_IJVector *v,
                                const int local_unknowns) {
  int start_ind = get_global_start_ind(local_unknowns);
  int end_ind = start_ind + local_unknowns - 1;
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, start_ind, end_ind, v);
  HYPRE_IJVectorSetObjectType(*v, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(*v);

  op_arg copy_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, DG_NP, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(v_dat->set, 1, copy_args);

  int setSize = v_dat->set->size;

  double *v_data;
  cudaMallocManaged(&v_data, DG_NP * setSize * sizeof(double));
  cudaMemcpy(v_data, v_dat->data_d, DG_NP * setSize * sizeof(double), cudaMemcpyDeviceToDevice);
  int *ind;
  cudaMallocManaged(&ind, DG_NP * setSize * sizeof(int));
  for(int i = 0; i < DG_NP * setSize; i++) {
    ind[i] = start_ind + i;
  }

  HYPRE_IJVectorSetValues(*v, DG_NP * setSize, ind, v_data);

  cudaFree(v_data);
  cudaFree(ind);
  op_mpi_set_dirtybit_cuda(1, copy_args);

  HYPRE_IJVectorAssemble(*v);
}

// Get a HYPRE vector as an OP2 dat
void HYPREUtils::vec_to_dat(op_dat v_dat, HYPRE_IJVector *v,
                            const int local_unknowns) {
  op_arg copy_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, DG_NP, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(v_dat->set, 1, copy_args);

  int setSize = v_dat->set->size;

  double *v_data;
  cudaMallocManaged(&v_data, DG_NP * setSize * sizeof(double));
  int *ind;
  cudaMallocManaged(&ind, DG_NP * setSize * sizeof(int));
  int start_ind = get_global_start_ind(local_unknowns);
  for(int i = 0; i < DG_NP * setSize; i++) {
    ind[i] = start_ind + i;
  }

  HYPRE_IJVectorGetValues(*v, DG_NP * setSize, ind, v_data);

  cudaMemcpy(v_dat->data_d, v_data, DG_NP * setSize * sizeof(double), cudaMemcpyDeviceToDevice);

  cudaFree(v_data);
  cudaFree(ind);
  op_mpi_set_dirtybit_cuda(1, copy_args);
}

void HYPREUtils::create_matrix(HYPRE_IJMatrix *mat, const int local_unknowns) {
  int start = 0;
  #ifdef INS_MPI
  start = get_global_mat_start_ind(local_unknowns);
  #endif
  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, start, start + local_unknowns - 1, start, start + local_unknowns - 1, mat);
  HYPRE_IJMatrixSetObjectType(*mat, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(*mat);
}
