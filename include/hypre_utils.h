#ifndef __INS_HYPRE_UTILS_H
#define __INS_HYPRE_UTILS_H

#include "op_seq.h"

#include "HYPRE_config.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

namespace HYPREUtils {
  void dat_to_new_vec(op_dat v_dat, HYPRE_IJVector *v, const int local_unknowns);
  void vec_to_dat(op_dat v_dat, HYPRE_IJVector *v, const int local_unknowns);

  void create_matrix(HYPRE_IJMatrix *mat, const int local_unknowns);
}

#endif
