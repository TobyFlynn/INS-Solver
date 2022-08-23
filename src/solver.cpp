#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

#include "solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_constants.h"
#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_operators.h"
#include "load_mesh.h"
#include "timing.h"

extern Timing *timer;
extern double reynolds, refVel;

using namespace std;

void Solver::reverse_vel() {
  int tmp = 1;
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&tmp, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[0][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[0][1], -1, OP_ID, DG_NP, "double", OP_WRITE));
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&tmp, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[1][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[1][1], -1, OP_ID, DG_NP, "double", OP_WRITE));
}

Solver::Solver(std::string filename, int prob) {
  problem = prob;

  // Ownership of the pointers is passed to DGMesh
  // so don't have to worry about freeing them
  double *coords_data;
  int *cells_data, *edge2node_data, *edge2cell_data, *bedge2node_data;
  int *bedge2cell_data, *bedge_type_data, *edgeNum_data, *bedgeNum_data;
  int numNodes_g, numCells_g, numEdges_g, numBoundaryEdges_g, numNodes;
  int numCells, numEdges, numBoundaryEdges;

  int pressure_dirichlet[3];
  int pressure_neumann[3];
  int viscosity_dirichlet[3];
  int viscosity_neumann[3];

  load_mesh(filename, &coords_data, &cells_data, &edge2node_data,
            &edge2cell_data, &bedge2node_data, &bedge2cell_data,
            &bedge_type_data, &edgeNum_data, &bedgeNum_data, &numNodes_g,
            &numCells_g, &numEdges_g, &numBoundaryEdges_g, &numNodes, &numCells,
            &numEdges, &numBoundaryEdges, pressure_dirichlet, pressure_neumann,
            viscosity_dirichlet, viscosity_neumann);

  mesh = new DGMesh(coords_data, cells_data, edge2node_data, edge2cell_data,
                    bedge2node_data, bedge2cell_data, bedge_type_data,
                    edgeNum_data, bedgeNum_data, numNodes_g, numCells_g,
                    numEdges_g, numBoundaryEdges_g, numNodes, numCells,
                    numEdges, numBoundaryEdges);

  data = new INSData(mesh);
  ls = new LS(mesh, data);
  pressurePoisson = new PetscPressureSolve(mesh, data, ls);
  viscosityPoisson = new PetscViscositySolve(mesh, data, ls);
  pMultigrid = new PMultigrid(mesh);

  pressurePoisson->setDirichletBCs(pressure_dirichlet);
  pressurePoisson->setNeumannBCs(pressure_neumann);
  viscosityPoisson->setDirichletBCs(viscosity_dirichlet);
  viscosityPoisson->setNeumannBCs(viscosity_neumann);
  pMultigrid->setDirichletBCs(pressure_dirichlet);
  pMultigrid->setNeumannBCs(pressure_neumann);

  op_partition("" STRINGIFY(OP2_PARTITIONER), "KWAY", mesh->cells, mesh->edge2cells, NULL);

  mesh->init();
  data->init();
  ls->init();
  pressurePoisson->init();
  viscosityPoisson->init();
  pMultigrid->init();

  // Set initial conditions
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[0][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[0][1], -1, OP_ID, DG_NP, "double", OP_WRITE));
  op_par_loop(set_ic, "set_ic", mesh->cells,
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->Q[1][0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->Q[1][1], -1, OP_ID, DG_NP, "double", OP_WRITE));

  dt = numeric_limits<double>::max();
  op_par_loop(calc_dt, "calc_dt", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_MIN));
  dt = dt / (DG_ORDER * DG_ORDER * refVel);
  op_printf("dt: %g\n", dt);
}

Solver::~Solver() {
  delete pMultigrid;
  delete viscosityPoisson;
  delete pressurePoisson;
  delete ls;
  delete data;
  delete mesh;
}

void Solver::set_sub_cycling(int sub_cycles) {
  sub_cycle = sub_cycles > 0;
  num_sub_cycles = sub_cycles > 0 ? sub_cycles : 1;
  macro_dt = dt * num_sub_cycles;
  op_printf("macro dt: %g\n", macro_dt);
  if(sub_cycle)
    op_printf("Number of sub cycles: %d\n", num_sub_cycles);
}

void Solver::advection_non_linear(op_dat u, op_dat v, op_dat Nx, op_dat Ny, double t) {
  // Tensor product of velocity with itself
  op_par_loop(advection_flux, "advection_flux", mesh->cells,
              op_arg_dat(u, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(v, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->F[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[3], -1, OP_ID, DG_NP, "double", OP_WRITE));

  div(mesh, data->F[0], data->F[1], Nx);
  div(mesh, data->F[2], data->F[3], Ny);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, u, 0.0, data->gQ[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, v, 0.0, data->gQ[1]);

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(data->flux[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(data->flux[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Exchange values on edges between elements
  op_par_loop(advection_faces, "advection_faces", mesh->edges,
              op_arg_dat(mesh->order,     -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(data->gQ[0],     -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[1],     -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->flux[0],   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->flux[1],   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  // Enforce BCs
  op_par_loop(advection_bc, "advection_bc", mesh->bedges,
              op_arg_dat(mesh->order,       0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->x,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->y,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[0],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[1],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->flux[0],   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->flux[1],   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[0], 1.0, Nx);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[1], 1.0, Ny);
}

void Solver::advection_non_linear(op_dat u0, op_dat v0, op_dat u1, op_dat v1, op_dat Nx, op_dat Ny, double t) {
  // Tensor product of the two velocities
  op_par_loop(advection_flux2, "advection_flux2", mesh->cells,
              op_arg_dat(u0, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(v0, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u1, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(v1, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->F[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[2], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->F[3], -1, OP_ID, DG_NP, "double", OP_WRITE));
  
  div(mesh, data->F[0], data->F[1], Nx);
  div(mesh, data->F[2], data->F[3], Ny);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, u0, 0.0, data->gQ[0]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, v0, 0.0, data->gQ[1]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, u1, 0.0, data->gQ[2]);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, v1, 0.0, data->gQ[3]);

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(data->flux[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(data->flux[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));
  
  op_par_loop(advection_faces2, "advection_faces2", mesh->edges,
              op_arg_dat(mesh->order,     -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(data->gQ[0],     -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[1],     -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[2],     -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[3],     -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->flux[0],   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->flux[1],   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));
  
  // Not sure with how the boundary flux is being calculated - maybe need to pass to times to kernel?
  op_par_loop(advection_bc2, "advection_bc2", mesh->bedges,
              op_arg_dat(mesh->order,       0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->x,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->y,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[0],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[1],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[2],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->gQ[3],     0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->flux[0],   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->flux[1],   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  // op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[0], 1.0, Nx);
  // op2_gemv(mesh, false, 1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[1], 1.0, Ny);
  // Mult by -1 for sub cycling
  op2_gemv(mesh, false, -1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[0], -1.0, Nx);
  op2_gemv(mesh, false, -1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, data->flux[1], -1.0, Ny);
}

void Solver::sub_cycle_velocity(op_dat u, op_dat v, op_dat u_l, op_dat v_l, double t, int num_cycles) {
  double time = t;

  op_par_loop(sub_cycle_init, "sub_cycle_init", mesh->cells,
                op_arg_dat(u,   -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(v,   -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(u_l, -1, OP_ID, DG_NP, "double", OP_WRITE),
                op_arg_dat(v_l, -1, OP_ID, DG_NP, "double", OP_WRITE));

  for(int cycle = 0; cycle < num_cycles; cycle++) {
    int x = -1;
    op_par_loop(sub_cycle_set_rkQ, "sub_cycle_set_rkQ", mesh->cells,
                op_arg_gbl(&x,  1, "int", OP_READ),
                op_arg_gbl(&dt, 1, "double", OP_READ),
                op_arg_dat(u_l,  -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(v_l,  -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[1][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[1][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rkQ[0],   -1, OP_ID, DG_NP, "double", OP_RW),
                op_arg_dat(data->rkQ[1],   -1, OP_ID, DG_NP, "double", OP_RW));

    for(int j = 0; j < 3; j++) {
      if(j == 0)
        time = t + cycle * dt;
      else if(j == 1)
        time = t + cycle * dt + dt;
      else
        time = t + cycle * dt + 0.5 * dt;
      advection_non_linear(u, v, data->rkQ[0], data->rkQ[1], data->rk[j][0], data->rk[j][1], time);

      if(j != 2) {
        op_par_loop(sub_cycle_set_rkQ, "sub_cycle_set_rkQ", mesh->cells,
                op_arg_gbl(&x,  1, "int", OP_READ),
                op_arg_gbl(&dt, 1, "double", OP_READ),
                op_arg_dat(u_l,  -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(v_l,  -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[1][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[1][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rkQ[0],   -1, OP_ID, DG_NP, "double", OP_RW),
                op_arg_dat(data->rkQ[1],   -1, OP_ID, DG_NP, "double", OP_RW));
      }
    }

    op_par_loop(sub_cycle_update_Q, "sub_cycle_update_Q", mesh->cells,
                op_arg_gbl(&dt, 1, "double", OP_READ),
                op_arg_dat(data->rk[0][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[0][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[1][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[1][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[2][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->rk[2][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(u_l, -1, OP_ID, DG_NP, "double", OP_RW),
                op_arg_dat(v_l, -1, OP_ID, DG_NP, "double", OP_RW));
  }
}

// Calculate Nonlinear Terms
void Solver::advection(int currentInd, double a0, double a1, double b0,
                       double b1, double g0, double t) {
  if(sub_cycle) {
    sub_cycle_velocity(data->Q[currentInd][0], data->Q[currentInd][1], data->Q_l[currentInd][0], data->Q_l[currentInd][1], t, num_sub_cycles);
    sub_cycle_velocity(data->Q[(currentInd + 1) % 2][0], data->Q[(currentInd + 1) % 2][1], data->Q_l[(currentInd + 1) % 2][0], data->Q_l[(currentInd + 1) % 2][1], t - macro_dt, num_sub_cycles * 2);

    // Calculate the intermediate velocity values
    op_par_loop(sub_cycle_advection_intermediate_vel, "sub_cycle_advection_intermediate_vel", mesh->cells,
                op_arg_gbl(&a0, 1, "double", OP_READ),
                op_arg_gbl(&a1, 1, "double", OP_READ),
                op_arg_dat(data->Q_l[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->Q_l[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->Q_l[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->Q_l[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->QT[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
                op_arg_dat(data->QT[1], -1, OP_ID, DG_NP, "double", OP_WRITE));
    
    // Needed for pressure BCs
    advection_non_linear(data->Q[currentInd][0], data->Q[currentInd][1], data->N[currentInd][0], data->N[currentInd][1], t);
  } else {
    // Calculate flux values
    advection_non_linear(data->Q[currentInd][0], data->Q[currentInd][1], data->N[currentInd][0], data->N[currentInd][1], t);

    // Calculate the intermediate velocity values
    op_par_loop(advection_intermediate_vel, "advection_intermediate_vel", mesh->cells,
                op_arg_gbl(&a0, 1, "double", OP_READ),
                op_arg_gbl(&a1, 1, "double", OP_READ),
                op_arg_gbl(&b0, 1, "double", OP_READ),
                op_arg_gbl(&b1, 1, "double", OP_READ),
                op_arg_gbl(&dt, 1, "double", OP_READ),
                op_arg_dat(data->Q[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->Q[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->Q[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->Q[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->N[currentInd][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->N[currentInd][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->N[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->N[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(data->QT[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
                op_arg_dat(data->QT[1], -1, OP_ID, DG_NP, "double", OP_WRITE));
  }
}

bool Solver::pressure(int currentInd, double a0, double a1, double b0,
                      double b1, double g0, double t) {
  timer->startTimer("Pressure Setup");

  div_with_central_flux(mesh, data->QT[0], data->QT[1], data->divVelT);
  curl(mesh, data->Q[currentInd][0], data->Q[currentInd][1], data->curlVel);
  grad(mesh, data->curlVel, data->gradCurlVel[0], data->gradCurlVel[1]);

  // Apply pressure boundary conditions
  op_par_loop(pressure_bc, "pressure_bc", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->order, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_gbl(&t, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->nx, 0, mesh->bedge2cells, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(mesh->ny, 0, mesh->bedge2cells, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(data->N[currentInd][0], 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->N[currentInd][1], 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->gradCurlVel[0], 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->gradCurlVel[1], 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->rho, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], 0, mesh->bedge2cells, 3 * DG_NPF, "double", OP_INC));

  if(problem == 1) {
    // TODO - update for p-adaptivity
    op_par_loop(pressure_bc2, "pressure_bc2", mesh->bedges,
                op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
                op_arg_gbl(&t, 1, "double", OP_READ),
                op_arg_gbl(&problem, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                op_arg_dat(data->prBC, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));
  }

  // Calculate RHS of pressure solve
  // This assumes that the boundaries will always be order DG_ORDER
  op_par_loop(pressure_rhs, "pressure_rhs", mesh->cells,
              op_arg_gbl(&b0, 1, "double", OP_READ),
              op_arg_gbl(&b1, 1, "double", OP_READ),
              op_arg_gbl(&macro_dt, 1, "double", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(data->dPdN[currentInd], -1, OP_ID, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, 3 * DG_NPF, "double", OP_RW),
              op_arg_dat(data->divVelT, -1, OP_ID, DG_NP, "double", OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, data->dPdN[(currentInd + 1) % 2], 1.0, data->divVelT);
  op2_gemv(mesh, false, 1.0, DGConstants::MASS, data->divVelT, 0.0, data->pRHS);
  timer->endTimer("Pressure Setup");

  // Call PETSc linear solver
  timer->startTimer("Pressure Linear Solve");
  pressurePoisson->setup();
  pressurePoisson->setBCValues(data->prBC);
  bool converged = pressurePoisson->solve(data->pRHS, data->p);

  // pMultigrid->setBCValues(data->prBC);
  // pMultigrid->set_rho(data->rho);
  // bool converged = pMultigrid->solve(data->pRHS, data->p);
  timer->endTimer("Pressure Linear Solve");

  // Calculate gradient of pressure
  grad_with_central_flux(mesh, data->p, data->dpdx, data->dpdy);

  // Calculate new velocity intermediate values
  op_par_loop(pressure_update_vel, "pressure_update_vel", mesh->cells,
              op_arg_gbl(&macro_dt, 1, "double", OP_READ),
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->dpdx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->dpdy, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QT[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QT[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[0], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->QTT[1], -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->dPdN[(currentInd + 1) % 2], -1, OP_ID, 3 * DG_NPF, "double", OP_WRITE),
              op_arg_dat(data->prBC, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  return converged;
}

bool Solver::viscosity(int currentInd, double a0, double a1, double b0,
                       double b1, double g0, double t) {
  timer->startTimer("Viscosity Setup");
  double time = t + macro_dt;

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(data->visBC[0], -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(data->visBC[1], -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Get BCs for viscosity solve
  op_par_loop(viscosity_bc, "viscosity_bc", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&time, 1, "double", OP_READ),
              op_arg_gbl(&problem, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->x, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->y, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(data->visBC[0], 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(data->visBC[1], 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  // Set up RHS for viscosity solve
  op_par_loop(viscosity_mm, "viscosity_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->visRHS[0], -1, OP_ID, DG_NP, "double", OP_WRITE));
  op_par_loop(viscosity_mm, "viscosity_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(data->QTT[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->visRHS[1], -1, OP_ID, DG_NP, "double", OP_WRITE));

  double factor = reynolds / macro_dt;
  op_par_loop(viscosity_rhs, "viscosity_rhs", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(data->rho,       -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->visRHS[0], -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(data->visRHS[1], -1, OP_ID, DG_NP, "double", OP_RW));

  timer->endTimer("Viscosity Setup");

  // Call PETSc linear solver
  timer->startTimer("Viscosity Linear Solve");
  factor = g0 * reynolds / macro_dt;
  viscosityPoisson->setup(factor);
  viscosityPoisson->setBCValues(data->visBC[0]);
  bool convergedX = viscosityPoisson->solve(data->visRHS[0], data->Q[(currentInd + 1) % 2][0]);

  viscosityPoisson->setBCValues(data->visBC[1]);
  bool convergedY = viscosityPoisson->solve(data->visRHS[1], data->Q[(currentInd + 1) % 2][1]);
  timer->endTimer("Viscosity Linear Solve");

  return convergedX && convergedY;
}

void Solver::update_surface(int currentInd) {
  timer->startTimer("Surface");
  ls->setVelField(data->Q[(currentInd + 1) % 2][0], data->Q[(currentInd + 1) % 2][1]);
  ls->step(dt, num_sub_cycles);
  timer->endTimer("Surface");
}

double Solver::getAvgPressureConvergance() {
  return pressurePoisson->getAverageConvergeIter();
}

double Solver::getAvgViscosityConvergance() {
  return viscosityPoisson->getAverageConvergeIter();
}
