#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

using namespace std;

Poisson_M::Poisson_M(DGMesh *m, INSData *d) : Poisson(m, d) {
  glb_ind_data   = (int *)calloc(mesh->numCells, sizeof(int));
  glb_indL_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indR_data  = (int *)calloc(mesh->numEdges, sizeof(int));
  glb_indBC_data = (int *)calloc(mesh->numBoundaryEdges, sizeof(int));
  op1_data       = (double *)calloc(15 * 15 * mesh->numCells, sizeof(double));
  op2_data[0]    = (double *)calloc(15 * 15 * mesh->numEdges, sizeof(double));
  op2_data[1]    = (double *)calloc(15 * 15 * mesh->numEdges, sizeof(double));
  op_bc_data     = (double *)calloc(7 * 15 * mesh->numBoundaryEdges, sizeof(double));

  glb_ind   = op_decl_dat(mesh->cells, 1, "int", glb_ind_data, "poisson_glb_ind");
  glb_indL  = op_decl_dat(mesh->edges, 1, "int", glb_indL_data, "poisson_glb_indL");
  glb_indR  = op_decl_dat(mesh->edges, 1, "int", glb_indR_data, "poisson_glb_indR");
  glb_indBC = op_decl_dat(mesh->bedges, 1, "int", glb_indBC_data, "poisson_glb_indBC");
  op1       = op_decl_dat(mesh->cells, 15 * 15, "double", op1_data, "poisson_op1");
  op2[0]    = op_decl_dat(mesh->edges, 15 * 15, "double", op2_data[0], "poisson_op20");
  op2[1]    = op_decl_dat(mesh->edges, 15 * 15, "double", op2_data[1], "poisson_op21");
  op_bc     = op_decl_dat(mesh->bedges, 7 * 15, "double", op_bc_data, "poisson_op_bc");
}

Poisson_M::~Poisson_M() {
  if(pMatInit)
    MatDestroy(&pMat);
  if(pMMatInit)
    MatDestroy(&pMMat);
  if(pBCMatInit)
    MatDestroy(&pBCMat);

  destroy_vec(&b);
  destroy_vec(&bc);
  destroy_vec(&rhs);
  destroy_vec(&x);

  MatDestroy(&op);
  KSPDestroy(&ksp);

  free(glb_ind_data);
  free(glb_indL_data);
  free(glb_indR_data);
  free(glb_indBC_data);
  free(op1_data);
  free(op2_data[0]);
  free(op2_data[1]);
}

void Poisson_M::init() {
  setGlbInd();
  op_par_loop(glb_ind_kernel, "glb_ind_kernel", mesh->edges,
              op_arg_dat(glb_ind, -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indL, -1, OP_ID, 1, "int", OP_WRITE),
              op_arg_dat(glb_indR, -1, OP_ID, 1, "int", OP_WRITE));
  op_par_loop(glb_ind_kernelBC, "glb_ind_kernelBC", mesh->bedges,
              op_arg_dat(glb_ind, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(glb_indBC, -1, OP_ID, 1, "int", OP_WRITE));
  setOp();
  setBCOP();
  createMatrix();
  createMassMatrix();
  createBCMatrix();

  create_vec(&b);
  create_vec(&bc, 21);
  create_vec(&rhs);
  create_vec(&x);

  create_mat(&op, 15 * mesh->cells->size, 15 * mesh->cells->size, 15 * 4);
  MatAssemblyBegin(op, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(op, MAT_FINAL_ASSEMBLY);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPCG);
  PC pc;
  KSPGetPC(ksp, &pc);
  // PCSetType(pc, PCICC);
  PCSetType(pc, PCNONE);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 1e4);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
}

bool Poisson_M::solve(op_dat b_dat, op_dat x_dat, bool addMass, double factor) {
  massMat = addMass;
  massFactor = factor;

  // Calculate RHS for linear solve by applying the BCs
  load_vec(&b, b_dat);
  load_vec(&bc, bc_dat, 21);
  MatMultAdd(pBCMat, bc, b, rhs);

  // Create matrix for linear solve, adding mass matrix scaled by a factor if required
  if(addMass) {
    MatCopy(pMat, op, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(op, factor, pMMat, DIFFERENT_NONZERO_PATTERN);
  } else {
    MatCopy(pMat, op, DIFFERENT_NONZERO_PATTERN);
  }

  KSPSetOperators(ksp, op, op);

  load_vec(&x, x_dat);

  // Solve
  KSPSolve(ksp, rhs, x);

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
    cout << "Number of iterations for linear solver: " << numIt << endl;
    cout << "Converged reason: " << reason << " Residual: " << residual << endl;
  }
  numberIter += numIt;
  solveCount++;

  // Get solution and free PETSc vectors and matrix
  Vec solution;
  KSPGetSolution(ksp, &solution);
  store_vec(&solution, x_dat);

  return converged;
}

void Poisson_M::setOp() {
  double tol = 1e-15;

  op_par_loop(poisson_mf2_op, "poisson_mf2_op", mesh->cells,
              op_arg_dat(data->cOP, -1, OP_ID, 15 * 15, "double", OP_READ),
              op_arg_gbl(&tol, 1, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, 15 * 15, "double", OP_WRITE));

  op_par_loop(poisson_mf2_opf, "poisson_mf2_opf", mesh->edges,
              op_arg_gbl(&tol, 1, "double", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->gOP[0], 0, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOP[1], 0, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOP[2], 0, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOPf[0], 0, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOPf[1], 0, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOPf[2], 0, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(op2[0], -1, OP_ID, 15 * 15, "double", OP_INC),
              op_arg_dat(op1, 0, mesh->edge2cells, 15 * 15, "double", OP_INC),
              op_arg_dat(data->gOP[0], 1, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOP[1], 1, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOP[2], 1, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOPf[0], 1, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOPf[1], 1, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOPf[2], 1, mesh->edge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, 15 * 15, "double", OP_INC),
              op_arg_dat(op1, 1, mesh->edge2cells, 15 * 15, "double", OP_INC));

  op_par_loop(poisson_mf2_opbf, "poisson_mf2_opbf", mesh->bedges,
              op_arg_gbl(&tol, 1, "double", OP_READ),
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(data->gOP[0], 0, mesh->bedge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOP[1], 0, mesh->bedge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(data->gOP[2], 0, mesh->bedge2cells, 15 * 15, "double", OP_READ),
              op_arg_dat(op1, 0, mesh->bedge2cells, 15 * 15, "double", OP_INC));
}

void Poisson_M::setBCOP() {
  double tol = 1e-15;
  // If not dirichlet BC, kernel will assume it is a neumann bc
  op_par_loop(poisson_mf2_bc, "poisson_mf2_bc", mesh->bedges,
              op_arg_gbl(&tol, 1, "double", OP_READ),
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[0], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[1], 1, "int", OP_READ),
              op_arg_gbl(&dirichlet[2], 1, "int", OP_READ),
              op_arg_dat(data->mD[0], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[1], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(data->mD[2], 0, mesh->bedge2cells, 7 * 15, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, 21, "double", OP_READ),
              op_arg_dat(data->tau, 0, mesh->bedge2cells, 3, "double", OP_READ),
              op_arg_dat(op_bc, -1, OP_ID, 7 * 15, "double", OP_INC));
}
