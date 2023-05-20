#include "linear_solvers/amgx_amg.h"

#include "matrices/poisson_coarse_matrix.h"
#include "timing.h"
#include "config.h"
#include "utils.h"

#include <iostream>
#include <type_traits>
#include <stdexcept>

extern Timing *timer;
extern Config *config;

AMGX_resources_handle amgx_res_handle;
AMGX_config_handle amgx_config_handle;

void print_callback(const char *msg, int length) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) { printf("%s", msg); }
}

AmgXAMGSolver::AmgXAMGSolver(DGMesh *m) {
  bc = nullptr;
  mesh = m;
  nullspace = false;

  AMGX_SAFE_CALL(AMGX_initialize());
  AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
  AMGX_SAFE_CALL(AMGX_install_signal_handler());
  AMGX_SAFE_CALL(AMGX_config_create_from_file(&amgx_config_handle, "/home/u1717021/Code/PhD/INS-Solver/config/config.amgx"));

  #ifdef INS_MPI
  amgx_comm = MPI_COMM_WORLD;
  int devices[] = {0};
  AMGX_resources_create(&amgx_res_handle, amgx_config_handle, &amgx_comm, 1, devices);
  #else
  AMGX_resources_create_simple(&amgx_res_handle, amgx_config_handle);
  #endif

  AMGX_vector_create(&rhs_amgx, amgx_res_handle, AMGX_mode_dDDI);
  AMGX_vector_create(&soln_amgx, amgx_res_handle, AMGX_mode_dDDI);
  AMGX_solver_create(&solver_amgx, amgx_res_handle, AMGX_mode_dDDI, amgx_config_handle);
}

AmgXAMGSolver::~AmgXAMGSolver() {
  // TODO destroy matrix
  AMGX_matrix_destroy(*amgx_mat);
  AMGX_solver_destroy(solver_amgx);
  AMGX_vector_destroy(rhs_amgx);
  AMGX_vector_destroy(soln_amgx);
  AMGX_finalize_plugins();
  AMGX_finalize();
}

bool AmgXAMGSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("AmgXAMGSolver - solve");
  if(dynamic_cast<PoissonCoarseMatrix*>(matrix) == nullptr) {
    throw std::runtime_error("AmgXAMGSolver matrix should be of type PoissonCoarseMatrix\n");
  }

  PoissonCoarseMatrix *coarse_mat = dynamic_cast<PoissonCoarseMatrix*>(matrix);

  coarse_mat->getAmgXMat(&amgx_mat);

  AMGX_vector_bind(rhs_amgx, *amgx_mat);
  AMGX_vector_bind(soln_amgx, *amgx_mat);

  DG_FP *rhs_ptr = getOP2PtrHost(rhs, OP_READ);
  DG_FP *ans_ptr = getOP2PtrHost(ans, OP_READ);

  AMGX_vector_upload(rhs_amgx,  coarse_mat->getUnknowns(), 1, rhs_ptr);
  AMGX_vector_upload(soln_amgx, coarse_mat->getUnknowns(), 1, ans_ptr);

  releaseOP2PtrHost(rhs, OP_READ, rhs_ptr);
  releaseOP2PtrHost(ans, OP_READ, ans_ptr);

  op_printf("Uploaded Vecs\n");

  AMGX_solver_setup(solver_amgx, *amgx_mat);

  op_printf("Set up solver\n");

  AMGX_solver_solve(solver_amgx, rhs_amgx, soln_amgx);

  op_printf("Solved\n");

  AMGX_SOLVE_STATUS status;
  AMGX_solver_get_status(solver_amgx, &status);
  if(status != AMGX_SOLVE_SUCCESS) {
    switch(status) {
      case AMGX_SOLVE_FAILED:
        throw std::runtime_error("AMGX solve failed");
      case AMGX_SOLVE_DIVERGED:
        throw std::runtime_error("AMGX solve diverged");
      default:
        throw std::runtime_error("AMGX solve failed, unrecognised status");
    }
  }

  DG_FP *ans_ptr_w = getOP2PtrHost(ans, OP_WRITE);
  AMGX_vector_download(soln_amgx, ans_ptr_w);
  releaseOP2PtrHost(ans, OP_WRITE, ans_ptr_w);

  timer->endTimer("AmgXAMGSolver - solve");

  return true;
}
