#include "hypre_utils.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"

int get_start_ind(const int local_unknowns) {
  return get_global_mat_start_ind(local_unknowns);
}
#else
int get_start_ind(const int local_unknowns) {
  return 0;
}
#endif

void HYPREUtils::create_vec(HYPRE_IJVector *v, const int local_unknowns) {
  int start = get_start_ind(local_unknowns);
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, start, start + local_unknowns - 1, v);
  HYPRE_IJVectorSetObjectType(*v, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(*v);
}

void HYPREUtils::create_matrix(HYPRE_IJMatrix *mat, const int local_unknowns) {
  int start = get_start_ind(local_unknowns);
  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, start, start + local_unknowns - 1, start, start + local_unknowns - 1, mat);
  HYPRE_IJMatrixSetObjectType(*mat, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(*mat);
}
