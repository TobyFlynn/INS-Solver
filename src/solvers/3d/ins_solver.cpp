#include "solvers/3d/ins_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"

#include "timing.h"
// #include "shocks.h"
#include "linear_solvers/petsc_amg.h"
#include "linear_solvers/petsc_block_jacobi.h"
#include "linear_solvers/petsc_pmultigrid.h"

#include <string>
#include <iostream>
#include <stdexcept>

extern Timing *timer;

INSSolver3D::INSSolver3D(DGMesh3D *m) {
  mesh = m;

  std::string name;
  DG_FP *dg_np_data = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  for(int i = 0; i < 3; i++) {
    name = "ins_solver_vel0" + std::to_string(i);
    vel[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_vel1" + std::to_string(i);
    vel[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_n0" + std::to_string(i);
    n[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_n1" + std::to_string(i);
    n[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_velT" + std::to_string(i);
    velT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
    name = "ins_solver_velTT" + std::to_string(i);
    velTT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
  }
  for(int i = 0; i < 9; i++) {
    name = "ins_solver_tmp_np" + std::to_string(i);
    tmp_np[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, name.c_str());
  }
  pr = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, dg_np_data, "ins_solver_pr");
  free(dg_np_data);

  DG_FP *dg_npf_data = (DG_FP *)calloc(4 * DG_NPF * mesh->cells->size, sizeof(DG_FP));
  tmp_npf[0] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_tmp_npf0");
  tmp_npf[1] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_tmp_npf1");
  tmp_npf[2] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_tmp_npf2");
  dPdN[0]    = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_dPdN0");
  dPdN[1]    = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, dg_npf_data, "ins_solver_dPdN1");
  free(dg_npf_data);

  DG_FP *dg_npf_bc_data = (DG_FP *)calloc(DG_NPF * mesh->bfaces->size, sizeof(DG_FP));
  tmp_npf_bc = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, dg_npf_bc_data, "ins_solver_tmp_npf_bc");
  free(dg_npf_bc_data);

  int *bc_1_data = (int *)calloc(mesh->bfaces->size, sizeof(int));
  tmp_bc_1 = op_decl_dat(mesh->bfaces, 1, "int", bc_1_data, "ins_solver_tmp_bc_1");
  bc_types = op_decl_dat(mesh->bfaces, 1, "int", bc_1_data, "ins_solver_bc_types");
  free(bc_1_data);

  DG_FP *cell_1_data = (DG_FP *)calloc(mesh->cells->size, sizeof(DG_FP));
  art_vis = op_decl_dat(mesh->cells, 1, DG_FP_STR, cell_1_data, "ins_solver_art_vis");
  free(cell_1_data);

  f[0][0] = tmp_np[0]; f[0][1] = tmp_np[1]; f[0][2] = tmp_np[2];
  f[1][0] = tmp_np[3]; f[1][1] = tmp_np[4]; f[1][2] = tmp_np[5];
  f[2][0] = tmp_np[6]; f[2][1] = tmp_np[7]; f[2][2] = tmp_np[8];

  divVelT     = tmp_np[0];
  curlVel[0]  = tmp_np[1];
  curlVel[1]  = tmp_np[2];
  curlVel[2]  = tmp_np[3];
  curl2Vel[0] = tmp_np[4];
  curl2Vel[1] = tmp_np[5];
  curl2Vel[2] = tmp_np[6];
  dpdx = tmp_np[0];
  dpdy = tmp_np[1];
  dpdz = tmp_np[2];
  vis_coeff = tmp_np[0];
  vis_mm = tmp_np[1];

  advec_flux[0] = tmp_npf[0];
  advec_flux[1] = tmp_npf[1];
  advec_flux[2] = tmp_npf[2];

  pr_bc_types  = tmp_bc_1;
  vis_bc_types = tmp_bc_1;
  pr_bc  = tmp_npf_bc;
  vis_bc = tmp_npf_bc;

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  dt = 0.0;
  time = 0.0;

  currentInd = 0;

  // pressureMatrix = new PoissonMatrix3D(mesh);
  pressureMatrix = new PoissonSemiMatrixFree3D(mesh);
  // viscosityMatrix = new MMPoissonMatrix3D(mesh);
  viscosityMatrix = new MMPoissonMatrixFree3D(mesh);
  // pressureSolver = new PETScAMGSolver(mesh);
  pressureSolver = new PETScPMultigrid(mesh);
  // pressureSolver = new PMultigridPoissonSolver(mesh);
  // viscositySolver = new PETScBlockJacobiSolver(mesh);
  // viscositySolver = new PETScAMGSolver(mesh);
  viscositySolver = new PETScInvMassSolver(mesh);

  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_nullspace(true);
  viscositySolver->set_matrix(viscosityMatrix);
  viscositySolver->set_nullspace(false);
}

INSSolver3D::~INSSolver3D() {
  delete pressureMatrix;
  delete viscosityMatrix;
  delete pressureSolver;
  delete viscositySolver;
}

void INSSolver3D::init(const DG_FP re, const DG_FP refVel) {
  timer->startTimer("INSSolver3D - Init");
  reynolds = re;

  // Set initial conditions
  op_par_loop(ins_3d_set_ic, "ins_3d_set_ic", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->z, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[1][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1) * refVel);
  // dt *= 1e-2;
  op_printf("INS dt is %g\n", dt);

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_bc_types, "ins_3d_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

  timer->endTimer("INSSolver3D - Init");
}

void INSSolver3D::step() {
  timer->startTimer("INSSolver3D - Advection");
  advection();
  timer->endTimer("INSSolver3D - Advection");

  timer->startTimer("INSSolver3D - Pressure");
  pressure();
  timer->endTimer("INSSolver3D - Pressure");

  // timer->startTimer("Shock Capturing");
  // shock_capturing();
  // timer->endTimer("Shock Capturing");

  timer->startTimer("INSSolver3D - Viscosity");
  viscosity();
  timer->endTimer("INSSolver3D - Viscosity");

  currentInd = (currentInd + 1) % 2;
  time += dt;
  g0 = 1.5;
  a0 = 2.0;
  a1 = -0.5;
  b0 = 2.0;
  b1 = -1.0;
}

void INSSolver3D::advection() {
  op_par_loop(ins_3d_advec_0, "ins_3d_advec_0", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0][0], f[0][1], f[0][2], n[currentInd][0]);
  mesh->div(f[1][0], f[1][1], f[1][2], n[currentInd][1]);
  mesh->div(f[2][0], f[2][1], f[2][2], n[currentInd][2]);

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(advec_flux[0], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_flux[1], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_flux[2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

  // Flux across faces
  op_par_loop(ins_3d_advec_1, "ins_3d_advec_1", mesh->faces,
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_flux[0], -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(advec_flux[1], -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(advec_flux[2], -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));

  // Boundary flux
  if(mesh->bface2cells) {
    op_par_loop(ins_3d_advec_2, "ins_3d_advec_2", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_flux[0], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(advec_flux[1], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(advec_flux[2], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, advec_flux[0], 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, advec_flux[1], 1.0, n[currentInd][1]);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, advec_flux[2], 1.0, n[currentInd][2]);

  // Calculate the intermediate velocity values
  op_par_loop(ins_3d_advec_3, "ins_3d_advec_3", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}

void INSSolver3D::pressure() {
  timer->startTimer("INSSolver3D - Pressure RHS");
  // mesh->div(velT[0], velT[1], velT[2], divVelT);
  mesh->div_with_central_flux(velT[0], velT[1], velT[2], divVelT);
  mesh->curl(vel[currentInd][0], vel[currentInd][1], vel[currentInd][2], curlVel[0], curlVel[1], curlVel[2]);
  mesh->curl(curlVel[0], curlVel[1], curlVel[2], curl2Vel[0], curl2Vel[1], curl2Vel[2]);

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_pr_0, "ins_3d_pr_0", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dPdN[currentInd], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));
  }

  op_par_loop(ins_3d_pr_1, "ins_3d_pr_1", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(divVelT, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT);
  mesh->mass(divVelT);

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_pr_2, "ins_3d_pr_2", mesh->bfaces,
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(pr_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(pr_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("INSSolver3D - Pressure RHS");

  timer->startTimer("INSSolver3D - Pressure Linear Solve");
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureMatrix->calc_mat();
  pressureSolver->set_bcs(pr_bc);
  bool converged = pressureSolver->solve(divVelT, pr);
  // if(!converged)
  //   throw std::runtime_error("\nPressure solve failed to converge\n");
  timer->endTimer("INSSolver3D - Pressure Linear Solve");

  timer->startTimer("INSSolver3D - Pressure Projection");
  mesh->grad_with_central_flux(pr, dpdx, dpdy, dpdz);

  op_par_loop(ins_3d_pr_3, "ins_3d_pr_3", mesh->cells,
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdz, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));
  timer->endTimer("INSSolver3D - Pressure Projection");
}

void INSSolver3D::viscosity() {
  timer->startTimer("INSSolver3D - Viscosity RHS");
  mesh->mass(velTT[0]);
  mesh->mass(velTT[1]);
  mesh->mass(velTT[2]);

  DG_FP factor = reynolds / dt;
  op_par_loop(ins_3d_vis_0, "ins_3d_vis_0", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  DG_FP vis_time = time + dt;
  if(mesh->bface2cells) {
    op_par_loop(ins_3d_vis_1, "ins_3d_vis_1", mesh->bfaces,
                op_arg_gbl(&vis_time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("INSSolver3D - Viscosity RHS");

  timer->startTimer("INSSolver3D - Viscosity Linear Solve");
  factor = g0 * reynolds / dt;

  if(factor != viscosityMatrix->get_factor()) {
    viscosityMatrix->set_factor(factor);
    viscosityMatrix->set_bc_types(vis_bc_types);
    // viscosityMatrix->calc_mat();
    viscositySolver->setFactor(1.0 / factor);
  }
  viscositySolver->set_bcs(vis_bc);
  bool convergedX = viscositySolver->solve(velTT[0], vel[(currentInd + 1) % 2][0]);
  if(!convergedX)
    throw std::runtime_error("\nViscosity X solve failed to converge\n");

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_vis_2, "ins_3d_vis_2", mesh->bfaces,
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  viscositySolver->set_bcs(vis_bc);
  bool convergedY = viscositySolver->solve(velTT[1], vel[(currentInd + 1) % 2][1]);
  if(!convergedY)
    throw std::runtime_error("\nViscosity Y solve failed to converge\n");

  viscositySolver->set_bcs(vis_bc);
  bool convergedZ = viscositySolver->solve(velTT[2], vel[(currentInd + 1) % 2][2]);
  if(!convergedZ)
    throw std::runtime_error("\nViscosity Z solve failed to converge\n");

  timer->endTimer("INSSolver3D - Viscosity Linear Solve");
}
/*
void INSSolver3D::shock_capturing() {
  // TODO maybe should be doing detector on something that isn't pressure?
  discontinuity_sensor(mesh, pr, art_vis, h);

  DG_FP factor = g0 * reynolds / dt;
  op_par_loop(ins_set_art_vis, "ins_set_art_vis", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(art_vis,   -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vis_coeff, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vis_mm,    -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  // const DG_FP *art_vis_ptr = (DG_FP *)art_vis->data;
  // std::cout << "Artificial Vis:" << std::endl;
  // for(int i = 0; i < mesh->cells->size; i++) {
  //   if(art_vis_ptr[i] > 0.0)
  //     std::cout << art_vis_ptr[i] << std::endl;
  // }
}
*/
DG_FP INSSolver3D::get_time() {
  return time;
}

DG_FP INSSolver3D::get_dt() {
  return dt;
}

void INSSolver3D::dump_data(const std::string &filename) {
  timer->startTimer("INSSolver3D - Dump Data");
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(mesh->z, filename.c_str());
  op_fetch_data_hdf5_file(vel[0][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][1], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][2], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][1], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][2], filename.c_str());
  op_fetch_data_hdf5_file(velT[0], filename.c_str());
  op_fetch_data_hdf5_file(velT[1], filename.c_str());
  op_fetch_data_hdf5_file(velT[2], filename.c_str());
  op_fetch_data_hdf5_file(velTT[0], filename.c_str());
  op_fetch_data_hdf5_file(velTT[1], filename.c_str());
  op_fetch_data_hdf5_file(velTT[2], filename.c_str());
  op_fetch_data_hdf5_file(pr, filename.c_str());
  timer->endTimer("INSSolver3D - Dump Data");
}
