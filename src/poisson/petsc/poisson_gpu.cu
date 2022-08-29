#include "petsc_poisson.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

#include "dg_utils.h"
#include "timing.h"

extern Timing *timer;

// Copy PETSc vec array to OP2 dat
void PetscPoissonSolve::copy_vec_to_dat(op_dat dat, const double *dat_d) {
  timer->startTimer("PETSc - vec2dat");
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
  timer->endTimer("PETSc - vec2dat");
}

// Copy OP2 dat to PETSc vec array
void PetscPoissonSolve::copy_dat_to_vec(op_dat dat, double *dat_d) {
  timer->startTimer("PETSc - dat2vec");
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
  timer->endTimer("PETSc - dat2vec");
}

// Create a PETSc vector for GPUs
void PetscPoissonSolve::create_vec(Vec *v) {
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECCUDA);
  VecSetSizes(*v, mat->unknowns, PETSC_DECIDE);
}

// Destroy a PETSc vector
void PetscPoissonSolve::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for GPUs
void PetscPoissonSolve::load_vec(Vec *v, op_dat v_dat) {
  double *v_ptr;
  VecCUDAGetArray(*v, &v_ptr);

  copy_dat_to_vec(v_dat, v_ptr);

  VecCUDARestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for GPUs
void PetscPoissonSolve::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecCUDAGetArrayRead(*v, &v_ptr);

  copy_vec_to_dat(v_dat, v_ptr);

  VecCUDARestoreArrayRead(*v, &v_ptr);
}

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  PetscPoissonSolve *poisson;
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

void PetscPoissonSolve::create_shell_mat() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreateShell(PETSC_COMM_WORLD, mat->unknowns, mat->unknowns, PETSC_DETERMINE, PETSC_DETERMINE, this, &pMat);
  MatShellSetOperation(pMat, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(pMat, VECCUDA);

  pMatInit = true;
}

PetscErrorCode precon(PC pc, Vec x, Vec y) {
  PetscPoissonSolve *poisson;
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

void PetscPoissonSolve::set_shell_pc(PC pc) {
  PCShellSetApply(pc, precon);
  PCShellSetContext(pc, this);
}

void PetscPoissonSolve::setMatrix() {
  timer->startTimer("PETSc - set mat");
  if(!pMatInit) {
    MatCreate(PETSC_COMM_WORLD, &pMat);
    pMatInit = true;
    MatSetSizes(pMat, mat->unknowns, mat->unknowns, PETSC_DECIDE, PETSC_DECIDE);

    #ifdef INS_MPI
    MatSetType(pMat, MATMPIAIJCUSPARSE);
    MatMPIAIJSetPreallocation(pMat, DG_NP * 4, NULL, 0, NULL);
    #else
    MatSetType(pMat, MATSEQAIJCUSPARSE);
    MatSeqAIJSetPreallocation(pMat, DG_NP * 4, NULL);
    #endif
    MatSetOption(pMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  }
  // Add cubature OP to Poisson matrix
  op_arg args[] = {
    op_arg_dat(mat->op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(mat->glb_ind, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->cells, 3, args);

  const int setSize = mesh->cells->size;
  double *op1_data = (double *)malloc(DG_NP * DG_NP * setSize * sizeof(double));
  int *glb   = (int *)malloc(setSize * sizeof(int));
  int *order = (int *)malloc(setSize * sizeof(int));
  cudaMemcpy(op1_data, mat->op1->data_d, setSize * DG_NP * DG_NP * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb, mat->glb_ind->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order, mesh->order->data_d, setSize * sizeof(int), cudaMemcpyDeviceToHost);
  op_mpi_set_dirtybit_cuda(3, args);

  MatSetOption(pMat, MAT_ROW_ORIENTED, PETSC_FALSE);

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
  }

  free(op1_data);
  free(glb);
  free(order);

  op_arg edge_args[] = {
    op_arg_dat(mat->op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(mat->op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
    op_arg_dat(mat->glb_indL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mat->glb_indR, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mat->orderL, -1, OP_ID, 1, "int", OP_READ),
    op_arg_dat(mat->orderR, -1, OP_ID, 1, "int", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(mesh->edges, 6, edge_args);
  double *op2L_data = (double *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(double));
  double *op2R_data = (double *)malloc(DG_NP * DG_NP * mesh->edges->size * sizeof(double));
  int *glb_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *glb_r = (int *)malloc(mesh->edges->size * sizeof(int));
  int *order_l = (int *)malloc(mesh->edges->size * sizeof(int));
  int *order_r = (int *)malloc(mesh->edges->size * sizeof(int));

  cudaMemcpy(op2L_data, mat->op2[0]->data_d, DG_NP * DG_NP * mesh->edges->size * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(op2R_data, mat->op2[1]->data_d, DG_NP * DG_NP * mesh->edges->size * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_l, mat->glb_indL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(glb_r, mat->glb_indR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_l, mat->orderL->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(order_r, mat->orderR->data_d, mesh->edges->size * sizeof(int), cudaMemcpyDeviceToHost);

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
  }

  free(op2L_data);
  free(op2R_data);
  free(glb_l);
  free(glb_r);
  free(order_l);
  free(order_r);

  op_mpi_set_dirtybit_cuda(6, edge_args);

  MatAssemblyBegin(pMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pMat, MAT_FINAL_ASSEMBLY);
  timer->endTimer("PETSc - set mat");
}
