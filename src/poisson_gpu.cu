#include "poisson.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

#include "dg_utils.h"

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

void PoissonSolve::create_shell_mat(Mat *m) {
  MatCreateShell(PETSC_COMM_WORLD, unknowns, unknowns, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(*m, VECCUDA);
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
  global_ind = get_global_start_index(glb_ind->set);
  #endif
  op_arg args[] = {
    op_arg_dat(glb_ind, -1, OP_ID, 1, "int", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 1, args);
  int *data_ptr = (int *)malloc(mesh->cells->size * sizeof(int));
  cudaMemcpy(data_ptr, glb_ind->data_d, glb_ind->set->size * sizeof(int), cudaMemcpyDeviceToHost);
  for(int i = 0; i < mesh->cells->size; i++) {
    data_ptr[i] = global_ind + i;
  }
  cudaMemcpy(glb_ind->data_d, data_ptr, glb_ind->set->size * sizeof(int), cudaMemcpyHostToDevice);
  op_mpi_set_dirtybit_cuda(1, args);
  free(data_ptr);
}

void PoissonSolve::setMatrix() {
  if(pMatInit) {
    MatDestroy(&pMat);
  }
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
  double *op1_data = (double *)malloc(DG_NP * DG_NP * mesh->cells->size * sizeof(double));
  int *glb   = (int *)malloc(mesh->cells->size * sizeof(int));
  int *order = (int *)malloc(mesh->cells->size * sizeof(int));
  cudaMemcpy(op1_data, op1->data_d, op1->set->size * DG_NP * DG_NP * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb, glb_ind->data_d, glb_ind->set->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order, mesh->order->data_d, mesh->order->set->size * sizeof(int), cudaMemcpyDeviceToHost);
  op_mpi_set_dirtybit_cuda(3, args);

  int currentRow = 0;
  int currentCol = 0;
  for(int i = 0; i < mesh->cells->size; i++) {
    int global_ind = glb[i];
    int Np, Nfp;
    DGUtils::basic_constants(order[i], &Np, &Nfp);

    // Convert data to row major format
    for(int m = 0; m < Np; m++) {
      for(int n = 0; n < Np; n++) {
        int row = currentRow + m;
        int col = currentCol + n;
        double val = op1_data[i * DG_NP * DG_NP + m + n * Np];
        MatSetValues(pMat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }

    glb[i] = currentRow;
    currentRow += Np;
    currentCol += Np;
  }

  free(op1_data);
  // free(glb);
  // free(order);

  op_arg edge_args[] = {
    op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->edges, 4, edge_args);
  double *op2L_data = (double *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(double));
  double *op2R_data = (double *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(double));
  int *glb_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *glb_r = (int *)malloc(mesh->edges->size * sizeof(int));

  cudaMemcpy(op2L_data, op2[0]->data_d, DG_NP * DG_NP * mesh->edges->size * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(op2R_data, op2[1]->data_d, DG_NP * DG_NP * mesh->edges->size * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_l, glb_indL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, glb_indR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);

  // Add Gauss OP and OPf to Poisson matrix
  for(int i = 0; i < mesh->edges->size; i++) {
    int leftElement = glb_l[i];
    int rightElement = glb_r[i];
    int leftRow  = glb[leftElement];
    int rightRow = glb[rightElement];
    int NpL, NpR, Nfp;
    DGUtils::basic_constants(order[leftElement], &NpL, &Nfp);
    DGUtils::basic_constants(order[rightElement], &NpR, &Nfp);

    // Gauss OPf
    // Convert data to row major format
    for(int m = 0; m < NpL; m++) {
      for(int n = 0; n < NpR; n++) {
        int row = leftRow + m;
        int col = rightRow + n;
        double val = op2L_data[i * DG_NP * DG_NP + m + n * NpL];
        MatSetValues(pMat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }
    // Convert data to row major format
    for(int m = 0; m < NpR; m++) {
      for(int n = 0; n < NpL; n++) {
        int row = rightRow + m;
        int col = leftRow + n;
        double val = op2R_data[i * DG_NP * DG_NP + m + n * NpR];
        MatSetValues(pMat, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
    }
  }

  free(op2L_data);
  free(op2R_data);
  free(glb_l);
  free(glb_r);

  free(glb);
  free(order);

  op_mpi_set_dirtybit_cuda(4, edge_args);

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);

}
