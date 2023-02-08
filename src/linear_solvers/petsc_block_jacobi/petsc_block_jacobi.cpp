#include "linear_solvers/petsc_block_jacobi.h"

#include "op_seq.h"

#include "utils.h"
#include "linear_solvers/petsc_utils.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

PETScBlockJacobiSolver::PETScBlockJacobiSolver(DGMesh2D *m) {
  nullspace = false;
  pMatInit = false;
  mesh = m;

  double *tmp_np = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  double *tmp_np_np = (double *)calloc(DG_NP * DG_NP * mesh->cells->size, sizeof(double));

  in  = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "block_jacobi_in");
  out = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, "block_jacobi_out");
  pre = op_decl_dat(mesh->cells, DG_NP * DG_NP, "double", tmp_np_np, "block_jacobi_pre");

  free(tmp_np_np);
  free(tmp_np);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 2.5e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCSHELL);
  set_shell_pc(pc);
  create_shell_mat();
  KSPSetOperators(ksp, pMat, pMat);
}

PETScBlockJacobiSolver::~PETScBlockJacobiSolver() {
  MatDestroy(&pMat);
  KSPDestroy(&ksp);
}

bool PETScBlockJacobiSolver::solve(op_dat rhs, op_dat ans) {
  // TODO only call when necessary
  calc_precond_mat();
  if(nullspace) {
    MatNullSpace ns;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &ns);
    MatSetNullSpace(pMat, ns);
    MatSetTransposeNullSpace(pMat, ns);
    MatNullSpaceDestroy(&ns);
  }

  matrix->apply_bc(rhs, bc);

  Vec b, x;
  PETScUtils::create_vec(&b, rhs->set);
  PETScUtils::create_vec(&x, ans->set);

  PETScUtils::load_vec(&b, rhs);
  PETScUtils::load_vec(&x, ans);

  KSPSolve(ksp, b, x);

  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  // Check that the solver converged
  bool converged = true;
  if(reason < 0) {
    double residual;
    KSPGetResidualNorm(ksp, &residual);
    converged = false;
    std::cout << "Number of iterations for linear solver: " << numIt << std::endl;
    std::cout << "Converged reason: " << reason << " Residual: " << residual << std::endl;
  }

  Vec solution;
  KSPGetSolution(ksp, &solution);
  PETScUtils::store_vec(&solution, ans);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  return converged;
}

void PETScBlockJacobiSolver::calc_rhs(const double *in_d, double *out_d) {
  // Copy u to OP2 dat
  PETScUtils::copy_vec_to_dat(in, in_d);

  matrix->mult(in, out);

  PETScUtils::copy_dat_to_vec(out, out_d);
}

// Matrix-free block-jacobi preconditioning function
void PETScBlockJacobiSolver::precond(const double *in_d, double *out_d) {
  PETScUtils::copy_vec_to_dat(in, in_d);

  op_par_loop(poisson_pre, "poisson_pre", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(pre, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_WRITE));

  PETScUtils::copy_dat_to_vec(out, out_d);
}

void PETScBlockJacobiSolver::calc_precond_mat() {
  const double *op1_ptr = getOP2PtrHost(matrix->op1, OP_READ);
  double *pre_ptr = getOP2PtrHost(pre, OP_WRITE);

  for(int i = 0; i < mesh->cells->size; i++) {
    const double *in_c = op1_ptr + i * matrix->op1->dim;
    double *inv_c      = pre_ptr + i * pre->dim;

    arma::mat a(in_c, DG_NP, DG_NP);
    arma::mat b(inv_c, DG_NP, DG_NP, false, true);

    b = arma::inv(a);
    // b = arma::inv_sympd(a);
  }

  releaseOP2PtrHost(matrix->op1, OP_READ, op1_ptr);
  releaseOP2PtrHost(pre, OP_WRITE, pre_ptr);
}
