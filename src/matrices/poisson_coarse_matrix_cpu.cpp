#include "matrices/poisson_coarse_matrix.h"

#include "op_seq.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

#include "utils.h"
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

#ifdef INS_MPI
#include "mpi.h"
#endif

#include "dg_mesh/dg_mesh_3d.h"

#ifdef INS_BUILD_WITH_HYPRE
void PoissonCoarseMatrix::setHYPREMatrix() {
  int global_size = getUnknowns();
  int local_size = getUnknowns();
  #ifdef INS_MPI
  global_size = global_sum(global_size);
  #endif
  // Keep track of how many non-zero entries locally
  // int nnz = 0;
  const int cell_set_size = _mesh->cells->size;
  const int faces_set_size = _mesh->faces->size + _mesh->faces->exec_size + _mesh->faces->nonexec_size;
  int nnz = cell_set_size * DG_NP_N1 * DG_NP_N1 + faces_set_size * DG_NP_N1 * DG_NP_N1 * 2;

  DG_FP *data_buf_ptr_h = (DG_FP *)malloc(nnz * sizeof(DG_FP));
  int *col_buf_ptr_h = (int *)malloc(nnz * sizeof(int));
  int *row_num_ptr_h = (int *)malloc(local_size * sizeof(int));
  int *num_col_ptr_h = (int *)malloc(local_size * sizeof(int));

  // Exchange halos
  DGMesh3D *mesh = dynamic_cast<DGMesh3D*>(_mesh);
  op_arg args[] = {
    op_arg_dat(op2[0], 0, mesh->flux2faces, op2[0]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[0], 1, mesh->flux2faces, op2[0]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[0], 2, mesh->flux2faces, op2[0]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[0], 3, mesh->flux2faces, op2[0]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[1], 0, mesh->flux2faces, op2[1]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[1], 1, mesh->flux2faces, op2[1]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[1], 2, mesh->flux2faces, op2[1]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[1], 3, mesh->flux2faces, op2[1]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(glb_indL, 0, mesh->flux2faces, glb_indL->dim, "int", OP_RW),
    op_arg_dat(glb_indL, 1, mesh->flux2faces, glb_indL->dim, "int", OP_RW),
    op_arg_dat(glb_indL, 2, mesh->flux2faces, glb_indL->dim, "int", OP_RW),
    op_arg_dat(glb_indL, 3, mesh->flux2faces, glb_indL->dim, "int", OP_RW),
    op_arg_dat(glb_indR, 0, mesh->flux2faces, glb_indR->dim, "int", OP_RW),
    op_arg_dat(glb_indR, 1, mesh->flux2faces, glb_indR->dim, "int", OP_RW),
    op_arg_dat(glb_indR, 2, mesh->flux2faces, glb_indR->dim, "int", OP_RW),
    op_arg_dat(glb_indR, 3, mesh->flux2faces, glb_indR->dim, "int", OP_RW)
  };
  op_mpi_halo_exchanges(mesh->fluxes, 16, args);
  op_mpi_wait_all(16, args);

  // Get data from OP2
  DG_FP *op1_data = getOP2PtrHost(op1, OP_READ);
  DG_FP *op2L_data = getOP2PtrHost(op2[0], OP_READ);
  DG_FP *op2R_data = getOP2PtrHost(op2[1], OP_READ);
  const int *glb   = (int *)glb_ind->data;
  const int *glb_l = (int *)glb_indL->data;
  const int *glb_r = (int *)glb_indR->data;

  if(!hypre_mat_init) {
    const int ilower = glb[0];
    const int iupper = glb[0] + local_size - 1;
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &hypre_mat);
    HYPRE_IJMatrixSetObjectType(hypre_mat, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(hypre_mat);
    hypre_mat_init = true;
  }

  std::map<int,std::vector<std::pair<int,DG_FP>>> mat_buffer;

  for(int c = 0; c < cell_set_size; c++) {
    // Add diagonal block to buffer
    int diag_base_col = glb[c];
    DG_FP *diag_data_ptr = op1_data + c * DG_NP_N1 * DG_NP_N1;
    for(int i = 0; i < DG_NP_N1; i++) {
      std::vector<std::pair<int,DG_FP>> row_buf;
      for(int j = 0; j < DG_NP_N1; j++) {
        int ind = i + j * DG_NP_N1;
        row_buf.push_back({diag_base_col + j, diag_data_ptr[ind]});
      }
      mat_buffer.insert({diag_base_col + i, row_buf});
    }
  }

  for(int k = 0; k < faces_set_size; k++) {
    if(glb_l[k] >= glb[0] && glb_l[k] < glb[0] + local_size) {
      int base_col = glb_r[k];
      DG_FP *face_data_ptr = op2L_data + k * DG_NP_N1 * DG_NP_N1;
      for(int i = 0; i < DG_NP_N1; i++) {
        std::vector<std::pair<int,DG_FP>> &row_buf = mat_buffer.at(glb_l[k] + i);
        for(int j = 0; j < DG_NP_N1; j++) {
          int ind = i + j * DG_NP_N1;
          row_buf.push_back({base_col + j, face_data_ptr[ind]});
        }
      }
    }
  }

  for(int k = 0; k < faces_set_size; k++) {
    if(glb_r[k] >= glb[0] && glb_r[k] < glb[0] + local_size) {
      int base_col = glb_l[k];
      DG_FP *face_data_ptr = op2R_data + k * DG_NP_N1 * DG_NP_N1;
      for(int i = 0; i < DG_NP_N1; i++) {
        std::vector<std::pair<int,DG_FP>> &row_buf = mat_buffer.at(glb_r[k] + i);
        for(int j = 0; j < DG_NP_N1; j++) {
          int ind = i + j * DG_NP_N1;
          row_buf.push_back({base_col + j, face_data_ptr[ind]});
        }
      }
    }
  }

  int current_nnz = 0;
  int current_row = 0;
  const int ilower = glb[0];
  const int iupper = glb[0] + local_size - 1;
  for(auto it = mat_buffer.begin(); it != mat_buffer.end(); it++) {
    std::sort(it->second.begin(), it->second.end());

    row_num_ptr_h[current_row] = it->first;
    int num_this_col = 0;
    for(int i = 0; i < it->second.size(); i++) {
      if(fabs(it->second[i].second) > 1e-8) {
        col_buf_ptr_h[current_nnz] = it->second[i].first;
        data_buf_ptr_h[current_nnz] = it->second[i].second;
        num_this_col++;
        current_nnz++;
      }
    }
    num_col_ptr_h[current_row] = num_this_col;
    current_row++;
  }


  HYPRE_IJMatrixSetValues(hypre_mat, local_size, num_col_ptr_h, row_num_ptr_h, col_buf_ptr_h, data_buf_ptr_h);

  free(data_buf_ptr_h);
  free(col_buf_ptr_h);
  free(row_num_ptr_h);
  free(num_col_ptr_h);

  releaseOP2PtrHost(op1, OP_READ, op1_data);
  releaseOP2PtrHost(op2[0], OP_READ, op2L_data);
  releaseOP2PtrHost(op2[1], OP_READ, op2R_data);

  HYPRE_IJMatrixAssemble(hypre_mat);
}
#endif
