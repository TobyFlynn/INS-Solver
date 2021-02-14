#include "poisson_rhs_solve.h"

#include <iostream>

#include "petscksp.h"
#include "petscvec.h"

#include "poisson_rhs.h"

using namespace std;

PetscErrorCode matAMult(Mat A, Vec x, Vec y);

void poisson_rhs_solve() {
  // Copy RHS to b Vec
  Vec b;
  VecCreateSeq(PETSC_COMM_SELF, 15 * data->numCells, &b);
  double *b_d;
  VecGetArray(b, &b_d);
  op_arg b_vec_petsc_args[] = {
    op_arg_dat(data->rhs, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges(data->cells, 1, b_vec_petsc_args);
  memcpy(b_d, (double *)data->rhs->data, 15 * data->numCells * sizeof(double));
  op_mpi_set_dirtybit(1, b_vec_petsc_args);
  VecRestoreArray(b, &b_d);

  Vec x;
  VecCreateSeq(PETSC_COMM_SELF, 15 * data->numCells, &x);

  // Create shell matrix
  Mat Amat;
  MatCreateShell(PETSC_COMM_SELF, 15 * data->numCells, 15 * data->numCells, PETSC_DETERMINE, PETSC_DETERMINE, NULL, &Amat);
  MatShellSetOperation(Amat, MATOP_MULT, (void(*)(void))matAMult);

  // Create KSP
  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetType(ksp, KSPGMRES);
  KSPSetOperators(ksp, Amat, Amat);

  // Solve
  KSPSolve(ksp, b, x);
  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  cout << "Number of iterations for linear solver: " << numIt << endl;
  cout << "Converged reason: " << reason << endl;

  // Get solution
  Vec solution;
  KSPGetSolution(ksp, &solution);
}

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  // Get PETSC data on GPU
  const double *x_d;
  double *y_d;
  VecGetArrayRead(x, &x_d);
  VecGetArray(y, &y_d);

  poisson_rhs(x_d, y_d);

  // Release PETSC data
  VecRestoreArrayRead(x, &x_d);
  VecRestoreArray(y, &y_d);

  return 0;
}
