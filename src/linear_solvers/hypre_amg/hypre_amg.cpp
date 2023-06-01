#include "linear_solvers/hypre_amg.h"

#include "matrices/poisson_coarse_matrix.h"
#include "timing.h"
#include "config.h"
#include "utils.h"

#include <iostream>
#include <type_traits>
#include <stdexcept>
#ifdef INS_CUDA
#include <cuda_runtime.h>
#endif

extern Timing *timer;
extern Config *config;

HYPREAMGSolver::HYPREAMGSolver(DGMesh *m) {
  bc = nullptr;
  mesh = m;
  nullspace = false;
  vec_init = false;

  HYPRE_Init();
  #if defined(HYPRE_USING_GPU)
  HYPRE_PrintDeviceInfo();
  HYPRE_SetSpGemmUseVendor(0);
  HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
  HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
  HYPRE_SetUseGpuRand(1);
  #endif
}

HYPREAMGSolver::~HYPREAMGSolver() {
  if(vec_init) {
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);
    /* Destroy solver and preconditioner */
    HYPRE_ParCSRPCGDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);
  }
  HYPRE_Finalize();
}

bool HYPREAMGSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("HYPREAMGSolver - solve");
  if(dynamic_cast<PoissonCoarseMatrix*>(matrix) == nullptr) {
    throw std::runtime_error("HYPREAMGSolver matrix should be of type PoissonCoarseMatrix\n");
  }

  PoissonCoarseMatrix *coarse_mat = dynamic_cast<PoissonCoarseMatrix*>(matrix);

  HYPRE_ParCSRMatrix *hypre_mat;
  bool setup_solver = false;
  timer->startTimer("HYPREAMGSolver - Get matrix and setup");
  if(coarse_mat->getHYPREMat(&hypre_mat)) {
    if(!vec_init) {
      int ilower, iupper, jlower, jupper;
      coarse_mat->getHYPRERanges(&ilower, &iupper, &jlower, &jupper);
      HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
      HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);

      HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
      HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);

      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
      HYPRE_PCGSetMaxIter(solver, 100);
      HYPRE_PCGSetTol(solver, 1e-7);
      HYPRE_PCGSetTwoNorm(solver, 1);
      HYPRE_PCGSetPrintLevel(solver, 2);
      HYPRE_PCGSetLogging(solver, 1);

      HYPRE_BoomerAMGCreate(&precond);
      HYPRE_BoomerAMGSetPrintLevel(precond, 1);
      HYPRE_BoomerAMGSetCoarsenType(precond, 8);
      // HYPRE_BoomerAMGSetOldDefault(precond);
      HYPRE_BoomerAMGSetRelaxType(precond, 18);
      HYPRE_BoomerAMGSetNumSweeps(precond, 1);
      HYPRE_BoomerAMGSetTol(precond, 0.0);
      HYPRE_BoomerAMGSetMaxIter(precond, 1);
      HYPRE_BoomerAMGSetKeepTranspose(precond, 1);
      HYPRE_BoomerAMGSetRAP2(precond, 1);
      HYPRE_BoomerAMGSetModuleRAP2(precond, 1);
      HYPRE_BoomerAMGSetStrongThreshold(precond, 0.5);
      HYPRE_BoomerAMGSetTruncFactor(precond, 0.2);

      HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                              (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

      vec_init = true;
      setup_solver = true;
    }
  }
  timer->endTimer("HYPREAMGSolver - Get matrix and setup");

  HYPRE_IJVectorInitialize(b);
  HYPRE_IJVectorInitialize(x);

  int ilower, iupper, jlower, jupper;
  coarse_mat->getHYPRERanges(&ilower, &iupper, &jlower, &jupper);

  timer->startTimer("HYPREAMGSolver - Transfer vec");
  DG_FP *rhs_ptr = getOP2PtrHost(rhs, OP_READ);
  DG_FP *ans_ptr = getOP2PtrHost(ans, OP_READ);

  const int num_unknowns_l = coarse_mat->getUnknowns();
  DG_FP *rhs_ptr_h = (DG_FP *)malloc(num_unknowns_l * sizeof(DG_FP));
  DG_FP *ans_ptr_h = (DG_FP *)malloc(num_unknowns_l * sizeof(DG_FP));
  int *ind_ptr_h = (int *)malloc(num_unknowns_l * sizeof(int));

  for(int i = 0; i < mesh->cells->size; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      rhs_ptr_h[i * DG_NP_N1 + j] = rhs_ptr[i * DG_NP + j];
      ans_ptr_h[i * DG_NP_N1 + j] = ans_ptr[i * DG_NP + j];
      ind_ptr_h[i * DG_NP_N1 + j] = ilower + i * DG_NP_N1 + j;
    }
  }

  releaseOP2PtrHost(rhs, OP_READ, rhs_ptr);
  releaseOP2PtrHost(ans, OP_READ, ans_ptr);

  #ifdef INS_CUDA
  DG_FP *data_rhs_ptr;
  cudaMalloc(&data_rhs_ptr, num_unknowns_l * sizeof(DG_FP));
  DG_FP *data_ans_ptr;
  cudaMalloc(&data_ans_ptr, num_unknowns_l * sizeof(DG_FP));
  int *data_ind_ptr;
  cudaMalloc(&data_ind_ptr, num_unknowns_l * sizeof(int));
  cudaMemcpy(data_rhs_ptr, rhs_ptr_h, num_unknowns_l * sizeof(DG_FP), cudaMemcpyHostToDevice);
  cudaMemcpy(data_ans_ptr, ans_ptr_h, num_unknowns_l * sizeof(DG_FP), cudaMemcpyHostToDevice);
  cudaMemcpy(data_ind_ptr, ind_ptr_h, num_unknowns_l * sizeof(int), cudaMemcpyHostToDevice);

  // Maybe need to transfer to device
  HYPRE_IJVectorSetValues(b, num_unknowns_l, data_ind_ptr, data_rhs_ptr);
  HYPRE_IJVectorSetValues(x, num_unknowns_l, data_ind_ptr, data_ans_ptr);

  cudaFree(data_rhs_ptr);
  #else
  HYPRE_IJVectorSetValues(b, num_unknowns_l, ind_ptr_h, rhs_ptr_h);
  HYPRE_IJVectorSetValues(x, num_unknowns_l, ind_ptr_h, ans_ptr_h);
  #endif

  free(rhs_ptr_h);

  HYPRE_IJVectorAssemble(b);
  HYPRE_IJVectorGetObject(b, (void **) &par_b);

  HYPRE_IJVectorAssemble(x);
  HYPRE_IJVectorGetObject(x, (void **) &par_x);
  timer->endTimer("HYPREAMGSolver - Transfer vec");

  if(setup_solver) {
    timer->startTimer("HYPREAMGSolver - Setup solver");
    HYPRE_ParCSRPCGSetup(solver, *hypre_mat, par_b, par_x);
    timer->endTimer("HYPREAMGSolver - Setup solver");
  }

  timer->startTimer("HYPREAMGSolver - Run solver");
  HYPRE_ParCSRPCGSolve(solver, *hypre_mat, par_b, par_x);
  timer->endTimer("HYPREAMGSolver - Run solver");

  int num_iterations;
  double final_res_norm;
  HYPRE_PCGGetNumIterations(solver, &num_iterations);
  HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

  op_printf("\n");
  op_printf("Iterations = %d\n", num_iterations);
  op_printf("Final Relative Residual Norm = %e\n", final_res_norm);
  op_printf("\n");

  timer->startTimer("HYPREAMGSolver - Transfer vec");
  DG_FP *ans_ptr_w = getOP2PtrHost(ans, OP_WRITE);

  #ifdef INS_CUDA
  HYPRE_IJVectorGetValues(x, num_unknowns_l, data_ind_ptr, data_ans_ptr);
  cudaMemcpy(ans_ptr_h, data_ans_ptr, num_unknowns_l * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaFree(data_ans_ptr);
  cudaFree(data_ind_ptr);
  #else
  HYPRE_IJVectorGetValues(x, num_unknowns_l, ind_ptr_h, ans_ptr_h);
  #endif

  for(int i = 0; i < mesh->cells->size; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      ans_ptr_w[i * DG_NP + j] = ans_ptr_h[i * DG_NP_N1 + j];
    }
  }

  releaseOP2PtrHost(ans, OP_WRITE, ans_ptr_w);

  free(ans_ptr_h);
  free(ind_ptr_h);
  timer->endTimer("HYPREAMGSolver - Transfer vec");

  return true;
}

void HYPREAMGSolver::set_tol(const DG_FP r_tol, const DG_FP a_tol) {
  op_printf("TODO set_tol for HYPREAMGSolver\n");
}
