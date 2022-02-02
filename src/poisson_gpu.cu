#include "poisson.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

#include "dg_utils.h"
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>

using namespace std;

int PoissonSolve::get_local_unknowns() {
  op_arg op2_args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->order->set, 1, op2_args);
  const int setSize = mesh->order->set->size;
  int *tempOrder = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(tempOrder, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  int local_unkowns = 0;
  for(int i = 0; i < setSize; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(tempOrder[i], &Np, &Nfp);
    local_unkowns += Np;
  }
  free(tempOrder);
  op_mpi_set_dirtybit_cuda(1, op2_args);
  return local_unkowns;
}

// Copy PETSc vec array to OP2 dat
void PoissonSolve::copy_vec_to_dat(op_dat dat, const double *dat_d) {
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_WRITE),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 2, copy_args);

  int setSize = dat->set->size;
  int *tempOrder = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(tempOrder, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  int vec_ind = 0;
  int block_start = 0;
  int block_count = 0;
  for(int i = 0; i < setSize; i++) {
    const int N = tempOrder[i];

    if(N == DG_ORDER) {
      if(block_count == 0) {
        block_start = i;
        block_count++;
        continue;
      } else {
        block_count++;
        continue;
      }
    } else {
      if(block_count != 0) {
        double *block_start_dat_c = (double *)dat->data_d + block_start * dat->dim;
        cudaMemcpy(block_start_dat_c, dat_d + vec_ind, block_count * DG_NP * sizeof(double), cudaMemcpyDeviceToDevice);
        vec_ind += DG_NP * block_count;
      }
      block_count = 0;
    }

    double *v_c = (double *)dat->data_d + i * dat->dim;
    int Np, Nfp;
    DGUtils::basic_constants(N, &Np, &Nfp);

    cudaMemcpy(v_c, dat_d + vec_ind, Np * sizeof(double), cudaMemcpyDeviceToDevice);
    vec_ind += Np;
  }

  if(block_count != 0) {
    double *block_start_dat_c = (double *)dat->data_d + block_start * dat->dim;
    cudaMemcpy(block_start_dat_c, dat_d + vec_ind, block_count * DG_NP * sizeof(double), cudaMemcpyDeviceToDevice);
    vec_ind += DG_NP * block_count;
  }
  free(tempOrder);
  op_mpi_set_dirtybit_cuda(2, copy_args);
}

// Copy OP2 dat to PETSc vec array
void PoissonSolve::copy_dat_to_vec(op_dat dat, double *dat_d) {
  op_arg copy_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 2, copy_args);

  int setSize = dat->set->size;
  int *tempOrder = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(tempOrder, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  int vec_ind = 0;
  int block_start = 0;
  int block_count = 0;
  for(int i = 0; i < setSize; i++) {
    const int N = tempOrder[i];

    if(N == DG_ORDER) {
      if(block_count == 0) {
        block_start = i;
        block_count++;
        continue;
      } else {
        block_count++;
        continue;
      }
    } else {
      if(block_count != 0) {
        const double *block_start_dat_c = (double *)dat->data_d + block_start * dat->dim;
        cudaMemcpy(dat_d + vec_ind, block_start_dat_c, block_count * DG_NP * sizeof(double), cudaMemcpyDeviceToDevice);
        vec_ind += DG_NP * block_count;
      }
      block_count = 0;
    }

    const double *v_c = (double *)dat->data_d + i * dat->dim;
    int Np, Nfp;
    DGUtils::basic_constants(N, &Np, &Nfp);

    cudaMemcpy(dat_d + vec_ind, v_c, Np * sizeof(double), cudaMemcpyDeviceToDevice);
    vec_ind += Np;
  }

  if(block_count != 0) {
    const double *block_start_dat_c = (double *)dat->data_d + block_start * dat->dim;
    cudaMemcpy(dat_d + vec_ind, block_start_dat_c, block_count * DG_NP * sizeof(double), cudaMemcpyDeviceToDevice);
    vec_ind += DG_NP * block_count;
  }
  free(tempOrder);
  op_mpi_set_dirtybit_cuda(2, copy_args);
}

// Create a PETSc vector for GPUs
void PoissonSolve::create_vec(Vec *v) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECCUDA);
  VecSetSizes(*v, unknowns, PETSC_DECIDE);
}

// Destroy a PETSc vector
void PoissonSolve::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for GPUs
void PoissonSolve::load_vec(Vec *v, op_dat v_dat) {
  double *v_ptr;
  VecCUDAGetArray(*v, &v_ptr);

  copy_dat_to_vec(v_dat, v_ptr);

  VecCUDARestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for GPUs
void PoissonSolve::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecCUDAGetArrayRead(*v, &v_ptr);

  copy_vec_to_dat(v_dat, v_ptr);

  VecCUDARestoreArrayRead(*v, &v_ptr);
}

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  PoissonSolve *poisson;
  MatShellGetContext(A, &poisson);
  const double *x_ptr;
  double *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);
  return 0;
}

void PoissonSolve::create_shell_mat() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreateShell(PETSC_COMM_WORLD, unknowns, unknowns, PETSC_DETERMINE, PETSC_DETERMINE, this, &pMat);
  MatShellSetOperation(pMat, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(pMat, VECCUDA);

  pMatInit = true;
}

PetscErrorCode precon(PC pc, Vec x, Vec y) {
  PoissonSolve *poisson;
  PCShellGetContext(pc, (void **)&poisson);
  const double *x_ptr;
  double *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  poisson->precond(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);
  return 0;
}

void PoissonSolve::set_shell_pc(PC pc) {
  PCShellSetApply(pc, precon);
  PCShellSetContext(pc, this);
}

void PoissonSolve::setGlbInd() {
  int global_ind = 0;
  #ifdef INS_MPI
  global_ind = get_global_mat_start_ind(unknowns);
  #endif
  op_arg args[] = {
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 2, args);

  const int setSize = mesh->cells->size;
  int *tempOrder = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(tempOrder, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  int *data_ptr = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(data_ptr, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);

  int ind = global_ind;
  for(int i = 0; i < mesh->cells->size; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(tempOrder[i], &Np, &Nfp);
    data_ptr[i] = ind;
    ind += Np;
  }

  cudaMemcpy(glb_ind->data_d, data_ptr, setSize * sizeof(int), cudaMemcpyHostToDevice);

  op_mpi_set_dirtybit_cuda(2, args);
  free(data_ptr);
  free(tempOrder);
}

void PoissonSolve::setMatrix() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreate(PETSC_COMM_WORLD, &pMat);
  pMatInit = true;
  MatSetSizes(pMat, unknowns, unknowns, PETSC_DECIDE, PETSC_DECIDE);

  #ifdef INS_MPI
  MatSetType(pMat, MATMPIAIJCUSPARSE);
  MatMPIAIJSetPreallocation(pMat, DG_NP * 4, NULL, 0, NULL);
  #else
  MatSetType(pMat, MATSEQAIJCUSPARSE);
  MatSeqAIJSetPreallocation(pMat, DG_NP * 4, NULL);
  #endif
  MatSetOption(pMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 3, args);

  const int setSize = mesh->cells->size;
  double *op1_data = (double *)malloc(DG_NP * DG_NP * setSize * sizeof(double));
  int *glb   = (int *)malloc(setSize * sizeof(int));
  int *order = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(op1_data, op1->data_d, setSize * DG_NP * DG_NP * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb, glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  op_mpi_set_dirtybit_cuda(3, args);

  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_FALSE);

  ofstream out("out.txt");
  out << "%%MatrixMarket matrix coordinate real general" << endl;
  out << to_string(mesh->cells->size * DG_NP) << " " << to_string(mesh->cells->size * DG_NP) << " " << to_string(mesh->cells->size * DG_NP * DG_NP + 2 * mesh->edges->size * DG_NP * DG_NP) << endl;

  for(int i = 0; i < setSize; i++) {
    int Np, Nfp;
    DGUtils::basic_constants(order[i], &Np, &Nfp);
    int currentRow = glb[i];
    int currentCol = glb[i];
    int idxm[DG_NP], idxn[DG_NP];
    for(int n = 0; n < DG_NP; n++) {
      idxm[n] = currentRow + n;
      idxn[n] = currentCol + n;
    }

    MatSetValues(pMat, Np, idxm, Np, idxn, &op1_data[i * DG_NP * DG_NP], INSERT_VALUES);
    for(int n = 0; n < DG_NP; n++) {
      for(int m = 0; m < DG_NP; m++) {
        out << to_string(idxm[n]) << " " << to_string(idxn[m]) << " " << to_string(op1_data[i * DG_NP * DG_NP + m * DG_NP + n]) << endl;
      }
    }
  }

  free(op1_data);
  free(glb);
  free(order);

  op_arg edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->edges, 6, edge_args);
  double *op2L_data = (double *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(double));
  double *op2R_data = (double *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(double));
  int *glb_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *glb_r = (int *)malloc(mesh->edges->size * sizeof(int));
  int *order_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *order_r = (int *)malloc(mesh->edges->size * sizeof(int));

  cudaMemcpy(op2L_data, op2[0]->data_d, DG_NP * DG_NP * mesh->edges->size * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(op2R_data, op2[1]->data_d, DG_NP * DG_NP * mesh->edges->size * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_l, glb_indL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_l, orderL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_r, orderR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < mesh->edges->size; i++) {
    int leftRow = glb_l[i];
    int rightRow = glb_r[i];
    int NpL, NpR, Nfp;
    DGUtils::basic_constants(order_l[i], &NpL, &Nfp);
    DGUtils::basic_constants(order_r[i], &NpR, &Nfp);

    int idxl[DG_NP], idxr[DG_NP];
    for(int n = 0; n < DG_NP; n++) {
      idxl[n] = leftRow + n;
      idxr[n] = rightRow + n;
    }

    MatSetValues(pMat, NpL, idxl, NpR, idxr, &op2L_data[i * DG_NP * DG_NP], INSERT_VALUES);
    MatSetValues(pMat, NpR, idxr, NpL, idxl, &op2R_data[i * DG_NP * DG_NP], INSERT_VALUES);

    for(int n = 0; n < DG_NP; n++) {
      for(int m = 0; m < DG_NP; m++) {
        out << to_string(idxl[n]) << " " << to_string(idxr[m]) << " " << to_string(op2L_data[i * DG_NP * DG_NP + m * DG_NP + n]) << endl;
      }
    }

    for(int n = 0; n < DG_NP; n++) {
      for(int m = 0; m < DG_NP; m++) {
        out << to_string(idxr[n]) << " " << to_string(idxl[m]) << " " << to_string(op2R_data[i * DG_NP * DG_NP + m * DG_NP + n]) << endl;
      }
    }
  }

  out.close();

  free(op2L_data);
  free(op2R_data);
  free(glb_l);
  free(glb_r);
  free(order_l);
  free(order_r);

  op_mpi_set_dirtybit_cuda(6, edge_args);

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);

}

void PressureSolve::setAMGXMat() {
  const int mat_size_in_blocks = mesh->cells->size * DG_NP;
  const int num_blocks = mesh->cells->size * DG_NP * DG_NP + 2 * mesh->edges->size * DG_NP * DG_NP;
  const int block_size = 1;

  op_arg op2_cell_args[] = {
    op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ)
  };

  op_arg op2_edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ)
  };

  double *diag_mat_data_d, *ordered_mat_data_d;
  cudaMalloc((void**)&diag_mat_data_d, mesh->cells->size * DG_NP * DG_NP * sizeof(double));
  cudaMalloc((void**)&ordered_mat_data_d, num_blocks * block_size * sizeof(double));
  int *glb_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *glb_r = (int *)malloc(mesh->edges->size * sizeof(int));

  op_mpi_halo_exchanges_cuda(mesh->cells, 1, op2_cell_args);
  cudaMemcpy(diag_mat_data_d, op1->data_d, mesh->cells->size * DG_NP * DG_NP * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, op2_cell_args);
  op_mpi_halo_exchanges_cuda(mesh->edges, 4, op2_edge_args);
  cudaMemcpy(glb_l, glb_indL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);

  const double *unordered_op1_d  = diag_mat_data_d;
  const double *unordered_op2L_d = op2[0]->data_d;
  const double *unordered_op2R_d = op2[1]->data_d;

  int *row_ptrs = (int *)calloc(mat_size_in_blocks + 1, sizeof(int));
  int *col_ptrs = (int *)calloc(num_blocks, sizeof(int));

  int current_col_ptrs_index = 0;
  int current_amgx_row = 0;
  for(int current_cell = 0; current_cell < mesh->cells->size; current_cell++) {
    vector<pair<int, const double *>> data_ptrs;

    data_ptrs.push_back(make_pair(current_cell * DG_NP, &unordered_op1_d[current_cell * DG_NP * DG_NP]));

    // Search op2L and op2R for blocks on this row
    for(int i = 0; i < mesh->edges->size; i++) {
      if(glb_l[i] / DG_NP == current_cell && glb_l[i] != glb_r[i]) {
        data_ptrs.push_back(make_pair(glb_r[i], &unordered_op2L_d[i * DG_NP * DG_NP]));
      }

      if(glb_r[i] / DG_NP == current_cell && glb_l[i] != glb_r[i]) {
        data_ptrs.push_back(make_pair(glb_l[i], &unordered_op2R_d[i * DG_NP * DG_NP]));
      }
    }

    sort(data_ptrs.begin(), data_ptrs.end());

    for(int i = 0; i < DG_NP; i++) {
      row_ptrs[current_amgx_row] = current_col_ptrs_index;
      for(int n = 0; n < data_ptrs.size(); n++) {
        cudaMemcpy(ordered_mat_data_d + current_col_ptrs_index, data_ptrs[n].second + i * DG_NP, DG_NP * sizeof(double), cudaMemcpyDeviceToDevice);
        for(int col = data_ptrs[n].first; col < data_ptrs[n].first + DG_NP; col++) {
          col_ptrs[current_col_ptrs_index] = col;
          current_col_ptrs_index++;
        }
      }
      current_amgx_row++;
    }
  }

  op_mpi_set_dirtybit_cuda(4, op2_edge_args);

  row_ptrs[mat_size_in_blocks] = num_blocks;

  op_printf("Number of blocks: %i\ncol_ptr_index: %i\n", num_blocks, current_col_ptrs_index);

  int *row_ptrs_d, *col_ptrs_d;
  cudaMalloc((void**)&row_ptrs_d, (mat_size_in_blocks + 1) * sizeof(int));
  cudaMalloc((void**)&col_ptrs_d, num_blocks * sizeof(int));

  cudaMemcpy(row_ptrs_d, row_ptrs, (mat_size_in_blocks + 1) * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(col_ptrs_d, col_ptrs, num_blocks * sizeof(int), cudaMemcpyHostToDevice);

  free(col_ptrs);
  free(row_ptrs);

  free(glb_l);
  free(glb_r);

  AMGX_SAFE_CALL(AMGX_matrix_upload_all(matrix, mat_size_in_blocks, num_blocks, block_size, block_size, row_ptrs_d, col_ptrs_d, (void *)ordered_mat_data_d, nullptr));

  cudaFree(col_ptrs_d);
  cudaFree(row_ptrs_d);
  cudaFree(ordered_mat_data_d);
  cudaFree(diag_mat_data_d);
}

void PressureSolve::uploadAMGXVec(AMGX_vector_handle *vec, op_dat dat) {
  op_arg op2_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 1, op2_args);

  AMGX_SAFE_CALL(AMGX_vector_upload(*vec, mesh->cells->size * DG_NP, 1, (void *)dat->data_d));

  op_mpi_set_dirtybit_cuda(1, op2_args);
}

void PressureSolve::downloadAMGXVec(AMGX_vector_handle *vec, op_dat dat) {
  op_arg op2_args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_NP, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 1, op2_args);

  AMGX_SAFE_CALL(AMGX_vector_download(*vec, (void *)dat->data_d));

  op_mpi_set_dirtybit_cuda(1, op2_args);
}
