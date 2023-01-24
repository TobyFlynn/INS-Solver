#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

#include "euler.h"

#include "op_seq.h"

#include <string>
#include "dg_constants.h"
#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_operators.h"

double gamma_e;

Euler::Euler(std::string &filename) {
  mesh = new DGMesh(filename);

  std::string name;
  double *tmp_np = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  double *tmp_g_np = (double *)calloc(DG_G_NP * mesh->cells->size, sizeof(double));
  for(int i = 0; i < 4; i++) {
    name = "Q" + std::to_string(i);
    Q[i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name = "F" + std::to_string(i);
    F[i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name = "G" + std::to_string(i);
    G[i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name = "rk_wQ" + std::to_string(i);
    rk_wQ[i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name = "rk_RHSQ0" + std::to_string(i);
    rk_RHSQ[0][i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name = "rk_RHSQ1" + std::to_string(i);
    rk_RHSQ[1][i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name = "rk_RHSQ2" + std::to_string(i);
    rk_RHSQ[2][i] = op_decl_dat(mesh->cells, DG_NP, "double", tmp_np, name.c_str());
    name = "gQ" + std::to_string(i);
    gQ[i] = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_g_np, name.c_str());
    name = "gRHSQ" + std::to_string(i);
    gRHSQ[i] = op_decl_dat(mesh->cells, DG_G_NP, "double", tmp_g_np, name.c_str());
  }
  free(tmp_g_np);
  free(tmp_np);

  gamma_e = 1.4;
  op_decl_const(1, "double", &gamma_e);
  op_decl_const(3 * 5, "int", DG_CONSTANTS);
  op_decl_const(3 * 3 * DG_NPF, "int", FMASK);
  op_decl_const(3 * DG_CUB_NP, "double", cubW_g);
  op_decl_const(3 * DG_GF_NP, "double", gaussW_g);

  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->edge2cells, NULL);

  mesh->init();
  init();
}

Euler::~Euler() {
  delete mesh;
}

void Euler::init() {
  op_par_loop(euler_ic, "euler_ic", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  // dt = 1e-4;
  dt = 0.00125;
}

void Euler::step() {
  rhs(Q, rk_RHSQ[0]);

  op_par_loop(euler_wQ_0, "euler_wQ_0", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_wQ[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  rhs(rk_wQ, rk_RHSQ[1]);

  op_par_loop(euler_wQ_1, "euler_wQ_1", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_wQ[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(rk_wQ[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  rhs(rk_wQ, rk_RHSQ[2]);

  op_par_loop(euler_wQ_2, "euler_wQ_2", mesh->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(Q[1], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(Q[2], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(Q[3], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(rk_RHSQ[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[0][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[1][3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[2][0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[2][1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[2][2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk_RHSQ[2][3], -1, OP_ID, DG_NP, "double", OP_READ));
}

void Euler::rhs(op_dat *wQ, op_dat *RHSQ) {
  op_par_loop(euler_flux, "euler_flux", mesh->cells,
              op_arg_dat(wQ[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(wQ[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(wQ[2], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(wQ[3], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(F[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(F[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(F[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(F[3], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(G[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(G[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(G[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(G[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  for(int i = 0; i < 4; i++) {
    div_weak(mesh, F[i], G[i], RHSQ[i]);
    op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, wQ[i], 0.0, gQ[i]);
  }

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(gRHSQ[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(gRHSQ[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(gRHSQ[2], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(gRHSQ[3], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  op_par_loop(euler_edges, "euler_edges", mesh->edges,
              op_arg_dat(mesh->order,     -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gQ[0],    -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gQ[1],    -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gQ[2],    -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gQ[3],    -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gRHSQ[0], -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(gRHSQ[1], -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(gRHSQ[2], -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(gRHSQ[3], -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  for(int i = 0; i < 4; i++) {
    op2_gemv(mesh, false, -1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, gRHSQ[i], 1.0, RHSQ[i]);
    // op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, gRHSQ[i], 1.0, RHSQ[i]);
    // inv_mass(mesh, RHSQ[i]);
  }
}

void Euler::dump_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(Q[0], filename.c_str());
  op_fetch_data_hdf5_file(Q[1], filename.c_str());
  op_fetch_data_hdf5_file(Q[2], filename.c_str());
  op_fetch_data_hdf5_file(Q[3], filename.c_str());
}
