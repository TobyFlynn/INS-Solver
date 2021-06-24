#include "pressure_solve.h"

#include "op_seq.h"

PressureSolve::PressureSolve(INSData *d, CubatureData *c, GaussData *g) {
  data = d;
  cData = c;
  gData = g;

  u_data   = (double *)calloc(15 * data->numCells, sizeof(double));
  rhs_data = (double *)calloc(15 * data->numCells, sizeof(double));

  u   = op_decl_dat(data->cells, 15, "double", u_data, "poisson_u");
  rhs = op_decl_dat(data->cells, 15, "double", rhs_data, "poisson_rhs");
}

PressureSolve::~PressureSolve() {
  free(u_data);
  free(rhs_data);

  destroy_vec(&b);
  destroy_vec(&x);
  KSPDestroy(&ksp);
  MatDestroy(&Amat);
}

void PressureSolve::setDirichletBCs(int *d) {
  dirichlet = d;
}

void PressureSolve::setNeumannBCs(int *n) {
  neumann = n;
}

void PressureSolve::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double PressureSolve::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}

void PressureSolve::init() {
  create_vec(&b);
  create_vec(&x);
  create_shell_mat(&Amat);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  KSPSetOperators(ksp, Amat, Amat);
  KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 1e5);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
}

bool PressureSolve::solve(op_dat b_dat, op_dat x_dat) {
  return false;
}

void PressureSolve::calc_rhs(const double *u_d, double *rhs_d) {
  // Copy u to OP2 dat
  copy_u(u_d);

  

  copy_rhs(rhs_d);
}
