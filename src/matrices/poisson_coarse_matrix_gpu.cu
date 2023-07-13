#include "matrices/poisson_coarse_matrix.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

#include <vector>
#include <algorithm>
#include <map>

#include "utils.h"
#include "dg_utils.h"
#include "dg_global_constants/dg_global_constants_2d.h"

#include "timing.h"
extern Timing *timer;

int PoissonCoarseMatrix::getUnknowns() {
  const int setSize_ = _mesh->order->set->size;
  int unknowns = setSize_ * DG_NP_N1;
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
  op_mpi_halo_exchanges_grouped(_mesh->cells, 1, args, 2, 0);
  op_mpi_wait_all_grouped(1, args, 2, 0);

  const int setSize = _mesh->cells->size;
  int *data_ptr = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(data_ptr, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  #pragma omp parallel for
  for(int i = 0; i < _mesh->cells->size; i++) {
    data_ptr[i] = global_ind + i * DG_NP_N1;
  }

  cudaMemcpy(glb_ind->data_d, data_ptr, setSize * sizeof(int), cudaMemcpyHostToDevice);

  op_mpi_set_dirtybit_cuda(1, args);
  free(data_ptr);
}

void PoissonCoarseMatrix::setPETScMatrix() {
  if(!petscMatInit) {
    timer->startTimer("setPETScMatrix - Create Matrix");
    MatCreate(PETSC_COMM_WORLD, &pMat);
    petscMatInit = true;
    int unknowns = getUnknowns();
    MatSetSizes(pMat, unknowns, unknowns, PETSC_DECIDE, PETSC_DECIDE);

    #ifdef INS_MPI
    MatSetType(pMat, MATMPIAIJCUSPARSE);
    // Guess for the off diagonal preallocation at the moment
    MatMPIAIJSetPreallocation(pMat, DG_NP_N1 * (DG_NUM_FACES + 1), NULL, DG_NP_N1 * 2, NULL);
    #else
    MatSetType(pMat, MATSEQAIJCUSPARSE);
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

  timer->startTimer("setPETScMatrix - OP2 op1");
  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(op1, -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_grouped(_mesh->cells, 2, args, 2, 0);
  op_mpi_wait_all_grouped(2, args, 2, 0);
  timer->endTimer("setPETScMatrix - OP2 op1");

  timer->startTimer("setPETScMatrix - Copy op1 to host");
  #ifdef DG_OP2_SOA
  const int setSize = getSetSizeFromOpArg(&args[0]);
  #else
  const int setSize = _mesh->cells->size;
  #endif
  DG_FP *op1_data = (DG_FP *)malloc(DG_NP_N1 * DG_NP_N1 * setSize * sizeof(DG_FP));
  int *glb   = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(op1_data, op1->data_d, setSize * DG_NP_N1 * DG_NP_N1 * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  timer->endTimer("setPETScMatrix - Copy op1 to host");
  op_mpi_set_dirtybit_cuda(2, args);

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

    #ifdef DG_OP2_SOA
    DG_FP tmp_op1[DG_NP_N1 * DG_NP_N1];
    for(int n = 0; n < DG_NP_N1 * DG_NP_N1; n++) {
      tmp_op1[n] = op1_data[i + n * setSize];
    }
    MatSetValues(pMat, DG_NP_N1, idxm, DG_NP_N1, idxn, tmp_op1, INSERT_VALUES);
    #else
    MatSetValues(pMat, DG_NP_N1, idxm, DG_NP_N1, idxn, &op1_data[i * DG_NP_N1 * DG_NP_N1], INSERT_VALUES);
    #endif
  }
  timer->endTimer("setPETScMatrix - Set values op1");

  free(op1_data);
  free(glb);

  timer->startTimer("setPETScMatrix - OP2 op2");
  op_arg edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, DG_NP_N1 * DG_NP_N1, DG_FP_STR, OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_grouped(_mesh->faces, 4, edge_args, 2, 0);
  op_mpi_wait_all_grouped(4, edge_args, 2, 0);
  timer->endTimer("setPETScMatrix - OP2 op2");

  timer->startTimer("setPETScMatrix - Copy op2 to host");
  #ifdef DG_OP2_SOA
  const int faces_set_size = getSetSizeFromOpArg(&edge_args[0]);
  #else
  const int faces_set_size = _mesh->faces->size;
  #endif
  DG_FP *op2L_data = (DG_FP *)malloc(DG_NP_N1 * DG_NP_N1 * faces_set_size * sizeof(DG_FP));
  DG_FP *op2R_data = (DG_FP *)malloc(DG_NP_N1 * DG_NP_N1 * faces_set_size * sizeof(DG_FP));
  int *glb_l = (int *)malloc(faces_set_size * sizeof(int));
  int *glb_r = (int *)malloc(faces_set_size * sizeof(int));

  cudaMemcpy(op2L_data, op2[0]->data_d, DG_NP_N1 * DG_NP_N1 * faces_set_size * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaMemcpy(op2R_data, op2[1]->data_d, DG_NP_N1 * DG_NP_N1 * faces_set_size * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_l, glb_indL->data_d, faces_set_size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, faces_set_size * sizeof(int), cudaMemcpyDeviceToHost);
  timer->endTimer("setPETScMatrix - Copy op2 to host");

  op_mpi_set_dirtybit_cuda(4, edge_args);

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

    #ifdef DG_OP2_SOA
    DG_FP tmp_op2_l[DG_NP_N1 * DG_NP_N1], tmp_op2_r[DG_NP_N1 * DG_NP_N1];
    for(int n = 0; n < DG_NP_N1 * DG_NP_N1; n++) {
      tmp_op2_l[n] = op2L_data[i + n * faces_set_size];
      tmp_op2_r[n] = op2R_data[i + n * faces_set_size];
    }
    MatSetValues(pMat, DG_NP_N1, idxl, DG_NP_N1, idxr, tmp_op2_l, INSERT_VALUES);
    MatSetValues(pMat, DG_NP_N1, idxr, DG_NP_N1, idxl, tmp_op2_r, INSERT_VALUES);
    #else
    MatSetValues(pMat, DG_NP_N1, idxl, DG_NP_N1, idxr, &op2L_data[i * DG_NP_N1 * DG_NP_N1], INSERT_VALUES);
    MatSetValues(pMat, DG_NP_N1, idxr, DG_NP_N1, idxl, &op2R_data[i * DG_NP_N1 * DG_NP_N1], INSERT_VALUES);
    #endif
  }
  timer->endTimer("setPETScMatrix - Set values op2");

  free(op2L_data);
  free(op2R_data);
  free(glb_l);
  free(glb_r);

  timer->startTimer("setPETScMatrix - Assembly");
  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
  timer->endTimer("setPETScMatrix - Assembly");
}

#ifdef INS_MPI
#include "mpi.h"
#endif

#include "dg_mesh/dg_mesh_3d.h"

#ifdef INS_BUILD_WITH_AMGX

extern AMGX_resources_handle amgx_res_handle;
extern AMGX_config_handle amgx_config_handle;

void PoissonCoarseMatrix::setAmgXMatrix() {
  if(!amgx_mat_init) {
    AMGX_matrix_create(&amgx_mat, amgx_res_handle, AMGX_mode_dDDI);
  }

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
  // Which entry are on which rows
  int *row_ptr = (int *)malloc((local_size + 1) * sizeof(int));
  #ifdef INS_MPI
  int64_t *col_inds = (int64_t *)malloc(nnz * sizeof(int64_t));
  #else
  int *col_inds = (int *)malloc(nnz * sizeof(int));
  #endif
  DG_FP *data_ptr = (DG_FP *)malloc(nnz * sizeof(DG_FP));

  // Exchange halos
  DGMesh3D *mesh = dynamic_cast<DGMesh3D*>(_mesh);
  op_arg args[] = {
    op_arg_dat(op2[0], 0, mesh->face2cells, op2[0]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[0], 1, mesh->face2cells, op2[0]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[1], 0, mesh->face2cells, op2[1]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(op2[1], 1, mesh->face2cells, op2[1]->dim, DG_FP_STR, OP_RW),
    op_arg_dat(glb_indL, 0, mesh->face2cells, glb_indL->dim, "int", OP_RW),
    op_arg_dat(glb_indL, 1, mesh->face2cells, glb_indL->dim, "int", OP_RW),
    op_arg_dat(glb_indR, 0, mesh->face2cells, glb_indR->dim, "int", OP_RW),
    op_arg_dat(glb_indR, 1, mesh->face2cells, glb_indR->dim, "int", OP_RW)
  };
  op_mpi_halo_exchanges_grouped(mesh->faces, 8, args, 2, 0);
  op_mpi_wait_all_grouped(8, args, 2, 0);
  cudaDeviceSynchronize();

  // Get data from OP2
  DG_FP *op1_data = getOP2PtrHost(op1, OP_READ);
  DG_FP *op2L_data = getOP2PtrHost(op2[0], OP_READ);
  DG_FP *op2R_data = getOP2PtrHost(op2[1], OP_READ);
  int *glb   = (int *)malloc(cell_set_size * sizeof(int));
  cudaMemcpy(glb, glb_ind->data_d, cell_set_size * sizeof(int), cudaMemcpyDeviceToHost);
  int *glb_l = (int *)malloc(faces_set_size * sizeof(int));
  int *glb_r = (int *)malloc(faces_set_size * sizeof(int));
  cudaMemcpy(glb_l, glb_indL->data_d, faces_set_size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, faces_set_size * sizeof(int), cudaMemcpyDeviceToHost);

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
  for(auto it = mat_buffer.begin(); it != mat_buffer.end(); it++) {
    std::sort(it->second.begin(), it->second.end());

    row_ptr[current_row] = current_nnz;
    for(int i = 0; i < it->second.size(); i++) {
      if(fabs(it->second[i].second) > 1e-8) {
        col_inds[current_nnz] = it->second[i].first;
        data_ptr[current_nnz] = it->second[i].second;
        current_nnz++;
      }
    }
    current_row++;
  }
  row_ptr[current_row] = current_nnz;
  // printf("cr: %d ls: %d nnz: %d cnnz: %d\n", current_row, local_size, nnz, current_nnz);
  // printf("gs: %d ls: %d\n", global_size, local_size);

  releaseOP2PtrHost(op1, OP_READ, op1_data);
  releaseOP2PtrHost(op2[0], OP_READ, op2L_data);
  releaseOP2PtrHost(op2[1], OP_READ, op2R_data);
  free(glb);
  free(glb_l);
  free(glb_r);

  AMGX_SAFE_CALL(AMGX_pin_memory(row_ptr, (local_size + 1) * sizeof(int)));
  AMGX_SAFE_CALL(AMGX_pin_memory(col_inds, nnz * sizeof(int)));
  AMGX_SAFE_CALL(AMGX_pin_memory(data_ptr, nnz * sizeof(DG_FP)));

  if(!amgx_mat_init) {
    #ifdef INS_MPI
    int nrings;
    AMGX_config_get_default_number_of_rings(amgx_config_handle, &nrings);
    int nranks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // If no partition vector is given, we assume a partitioning with contiguous blocks (see example above). It is sufficient (and faster/more scalable)
    // to calculate the partition offsets and pass those into the API call instead of creating a full partition vector.
    int64_t* partition_offsets = (int64_t*)malloc((nranks+1) * sizeof(int64_t));
    // gather the number of rows on each rank, and perform an exclusive scan to get the offsets.
    int64_t n64 = local_size;
    partition_offsets[0] = 0; // rows of rank 0 always start at index 0
    MPI_Allgather(&n64, 1, MPI_INT64_T, &partition_offsets[1], 1, MPI_INT64_T, MPI_COMM_WORLD);
    for (int i = 2; i < nranks + 1; ++i) {
      partition_offsets[i] += partition_offsets[i-1];
    }
    global_size = partition_offsets[nranks]; // last element always has global number of rows

    // for(int i = 0; i < nranks + 1; i++) {
    //   op_printf("%d\n", partition_offsets[i]);
    // }

    AMGX_distribution_handle dist;
    AMGX_distribution_create(&dist, amgx_config_handle);
    AMGX_distribution_set_partition_data(dist, AMGX_DIST_PARTITION_OFFSETS, partition_offsets);
    AMGX_matrix_upload_distributed(amgx_mat, global_size, local_size, current_nnz, 1, 1, row_ptr, col_inds, data_ptr, NULL, dist);
    AMGX_distribution_destroy(dist);
    free(partition_offsets);
    // int nrings;
    // AMGX_config_get_default_number_of_rings(amgx_config_handle, &nrings);
    // AMGX_SAFE_CALL(AMGX_matrix_upload_all_global(amgx_mat, global_size, local_size, current_nnz, 1, 1,
    //                               row_ptr, col_inds, data_ptr, NULL, nrings, nrings, part_vec));
    #else
    AMGX_SAFE_CALL(AMGX_matrix_upload_all(amgx_mat, local_size, current_nnz, 1, 1, row_ptr, col_inds, data_ptr, NULL));
    #endif
    amgx_mat_init = true;
  } else {
    AMGX_matrix_replace_coefficients(amgx_mat, local_size, current_nnz, data_ptr, NULL);
  }

  AMGX_SAFE_CALL(AMGX_unpin_memory(row_ptr));
  AMGX_SAFE_CALL(AMGX_unpin_memory(col_inds));
  AMGX_SAFE_CALL(AMGX_unpin_memory(data_ptr));

  free(row_ptr);
  free(col_inds);
  free(data_ptr);
}
#endif

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
  op_arg args[] = {
    op_arg_dat(glb_indL, -1, OP_ID, glb_indL->dim, "int", OP_RW),
    op_arg_dat(glb_indR, -1, OP_ID, glb_indR->dim, "int", OP_RW)
  };
  op_mpi_halo_exchanges_grouped(_mesh->faces, 2, args, 2, 1);
  op_mpi_wait_all_grouped(2, args, 2, 1);
  cudaDeviceSynchronize();

  // Get data from OP2
  DG_FP *op1_data = getOP2PtrHostHE(op1, OP_READ);
  DG_FP *op2L_data = getOP2PtrHostHE(op2[0], OP_READ);
  DG_FP *op2R_data = getOP2PtrHostHE(op2[1], OP_READ);
  int *glb   = (int *)malloc(cell_set_size * sizeof(int));
  cudaMemcpy(glb, glb_ind->data_d, cell_set_size * sizeof(int), cudaMemcpyDeviceToHost);
  int *glb_l = (int *)malloc(faces_set_size * sizeof(int));
  int *glb_r = (int *)malloc(faces_set_size * sizeof(int));
  cudaMemcpy(glb_l, glb_indL->data_d, faces_set_size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, faces_set_size * sizeof(int), cudaMemcpyDeviceToHost);

  if(!hypre_mat_init) {
    const int ilower = glb[0];
    const int iupper = glb[0] + local_size - 1;
    // printf("ilower: %d iupper: %d\n", ilower, iupper);
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

  DG_FP *data_buf_ptr;
  cudaMalloc(&data_buf_ptr, current_nnz * sizeof(DG_FP));
  int *col_buf_ptr;
  cudaMalloc(&col_buf_ptr, current_nnz * sizeof(int));
  int *row_num_ptr;
  cudaMalloc(&row_num_ptr, local_size * sizeof(int));
  int *num_col_ptr;
  cudaMalloc(&num_col_ptr, local_size * sizeof(int));

  cudaMemcpy(row_num_ptr, row_num_ptr_h, local_size * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(num_col_ptr, num_col_ptr_h, local_size * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(col_buf_ptr, col_buf_ptr_h, current_nnz * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(data_buf_ptr, data_buf_ptr_h, current_nnz * sizeof(DG_FP), cudaMemcpyHostToDevice);

  HYPRE_IJMatrixSetValues(hypre_mat, local_size, num_col_ptr, row_num_ptr, col_buf_ptr, data_buf_ptr);

  cudaFree(data_buf_ptr);
  cudaFree(col_buf_ptr);
  cudaFree(row_num_ptr);
  cudaFree(num_col_ptr);
  free(data_buf_ptr_h);
  free(col_buf_ptr_h);
  free(row_num_ptr_h);
  free(num_col_ptr_h);

  releaseOP2PtrHostHE(op1, OP_READ, op1_data);
  releaseOP2PtrHostHE(op2[0], OP_READ, op2L_data);
  releaseOP2PtrHostHE(op2[1], OP_READ, op2R_data);
  free(glb);
  free(glb_l);
  free(glb_r);

  HYPRE_IJMatrixAssemble(hypre_mat);
  // HYPRE_IJMatrixPrint(hypre_mat, "IJ.out.hypre_mat");
  // op_printf("Finish HYPRE Assembly\n");
}
#endif
