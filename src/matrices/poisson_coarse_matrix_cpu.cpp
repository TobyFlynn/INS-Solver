#include "matrices/poisson_coarse_matrix.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

#include "dg_utils.h"
#include "dg_global_constants/dg_global_constants_2d.h"

#include "timing.h"
extern Timing *timer;

int PoissonCoarseMatrix::getUnknowns() {
  const int setSize = _mesh->order->set->size;
  int unknowns = setSize * DG_NP_N1;
  return unknowns;
}

void PoissonCoarseMatrix::set_glb_ind() {
  int unknowns = getUnknowns();
  int global_ind = 0;
  #ifdef INS_MPI
  global_ind = get_global_mat_start_ind(unknowns);
  #endif
  op_arg args[] = {
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges(_mesh->cells, 1, args);

  int *data_ptr = (int *)glb_ind->data;
  #pragma omp parallel for
  for(int i = 0; i < _mesh->cells->size; i++) {
    data_ptr[i] = global_ind + i * DG_NP_N1;
  }
  op_mpi_set_dirtybit(1, args);
}

void PoissonCoarseMatrix::setPETScMatrix() {
  if(!petscMatInit) {
    timer->startTimer("setPETScMatrix - Create Matrix");
    MatCreate(PETSC_COMM_WORLD, &pMat);
    petscMatInit = true;
    int unknowns = getUnknowns();
    MatSetSizes(pMat, unknowns, unknowns, PETSC_DECIDE, PETSC_DECIDE);

    #ifdef INS_MPI
    MatSetType(pMat, MATMPIAIJ);
    MatMPIAIJSetPreallocation(pMat, DG_NP_N1 * (DG_NUM_FACES + 1), NULL, DG_NP_N1 * 2, NULL);
    #else
    MatSetType(pMat, MATSEQAIJ);
    MatSeqAIJSetPreallocation(pMat, DG_NP_N1 * (DG_NUM_FACES + 1), NULL);
    #endif
    MatSetOption(pMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetOption(pMat, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);
    MatSetOption(pMat, MAT_STRUCTURAL_SYMMETRY_ETERNAL, PETSC_TRUE);
    MatSetOption(pMat, MAT_SPD, PETSC_TRUE);
    MatSetOption(pMat, MAT_SPD_ETERNAL, PETSC_TRUE);
    timer->endTimer("setPETScMatrix - Create Matrix");
  } else {
    MatSetOption(pMat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
    MatSetOption(pMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    // MatSetOption(pMat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
    MatSetOption(pMat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  }

  // Add cubature OP to Poisson matrix
  timer->startTimer("setPETScMatrix - OP2 op1");
  op_arg args[] = {
    op_arg_dat(op1, -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(_mesh->cells, 2, args);
  op_mpi_set_dirtybit(2, args);
  timer->endTimer("setPETScMatrix - OP2 op1");

  const DG_FP *op1_data = (DG_FP *)op1->data;
  const int *glb = (int *)glb_ind->data;

  #ifdef DG_COL_MAJ
  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_FALSE);
  #else
  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_TRUE);
  #endif

  timer->startTimer("setPETScMatrix - Set values op1");
  for(int i = 0; i < _mesh->cells->size; i++) {
    int currentRow = glb[i];
    int currentCol = glb[i];

    int idxm[DG_NP_N1], idxn[DG_NP_N1];
    for(int n = 0; n < DG_NP_N1; n++) {
      idxm[n] = currentRow + n;
      idxn[n] = currentCol + n;
    }

    MatSetValues(pMat, DG_NP_N1, idxm, DG_NP_N1, idxn, &op1_data[i * DG_NP_N1 * DG_NP_N1], INSERT_VALUES);
  }
  timer->endTimer("setPETScMatrix - Set values op1");

  timer->startTimer("setPETScMatrix - OP2 op2");
  op_arg edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges(_mesh->faces, 4, edge_args);
  op_mpi_set_dirtybit(4, edge_args);
  timer->endTimer("setPETScMatrix - OP2 op2");

  const DG_FP *op2L_data = (DG_FP *)op2[0]->data;
  const DG_FP *op2R_data = (DG_FP *)op2[1]->data;
  const int *glb_l = (int *)glb_indL->data;
  const int *glb_r = (int *)glb_indR->data;

  // Add Gauss OP and OPf to Poisson matrix
  timer->startTimer("setPETScMatrix - Set values op2");
  for(int i = 0; i < _mesh->faces->size; i++) {
    int leftRow = glb_l[i];
    int rightRow = glb_r[i];

    int idxl[DG_NP_N1], idxr[DG_NP_N1];
    for(int n = 0; n < DG_NP_N1; n++) {
      idxl[n] = leftRow + n;
      idxr[n] = rightRow + n;
    }

    MatSetValues(pMat, DG_NP_N1, idxl, DG_NP_N1, idxr, &op2L_data[i * DG_NP_N1 * DG_NP_N1], INSERT_VALUES);
    MatSetValues(pMat, DG_NP_N1, idxr, DG_NP_N1, idxl, &op2R_data[i * DG_NP_N1 * DG_NP_N1], INSERT_VALUES);
  }
  timer->endTimer("setPETScMatrix - Set values op2");

  timer->startTimer("setPETScMatrix - Assembly");
  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
  timer->endTimer("setPETScMatrix - Assembly");
}
