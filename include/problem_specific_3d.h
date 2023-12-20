#ifndef __INS_PROBLEM_SPECIFIC_3D_H
#define __INS_PROBLEM_SPECIFIC_3D_H

#define EXAMPLE_3D_TGV

#if defined(EXAMPLE_3D_TGV)
#include "examples/3d_tgv.h"
#elif defined(EXAMPLE_3D_LIQUID_WHISTLE)
#include "examples/3d_liquid_whistle.h"
#else

#if defined(INS_CUDA)
#define DEVICE_PREFIX __device__
#elif defined(INS_HIP)
#define DEVICE_PREFIX __device__
#else
#define DEVICE_PREFIX inline
#endif

// BC types for linear solvers
#define BC_DIRICHLET 0
#define BC_NEUMANN 1

// Hardcoded BC types, do not edit
#define BC_TYPE_NO_SLIP 0
#define BC_TYPE_SLIP 1
#define BC_TYPE_NATURAL_OUTFLOW 2
// Add custom BC types below (number must be greater than 0), for example:
#define BC_TYPE_INFLOW 3

// Required definitions
#define LS_CAP 1.0
// Problem specifc definitions
#define LW_LENGTH 50.0
#define LW_LENGTH_SHORT 10.0
#define LW_INLET_RADIUS 0.5
#define LW_INLET_NO_SLIP_RADIUS 0.75
#define LW_INLET_LENGTH 2.0
#define LW_RADIUS 5.833
#define LW_BLADE_START_X 0.6

/************************************************************************
 * You can edit the body of the functions below but not their signature *
 ************************************************************************/

// Set the initial conditions of the problem
DEVICE_PREFIX void ps3d_set_ic(const DG_FP x, const DG_FP y, const DG_FP z,
                               DG_FP &u, DG_FP &v, DG_FP &w) {
  // Liquid Whistle Initial Conditions
  // DG_FP tmp = y * y + z * z - LW_INLET_RADIUS * LW_INLET_RADIUS;
  // tmp = fabs(tmp);
  // tmp = fmin(1.0, tmp / 0.1);
  // u = tmp * fmax(0.0, 0.0 - (x + 1.0));
  // v = 0.0;
  // w = 0.0;

  // TGV Initial Conditions
  u = sin(x) * cos(y) * cos(z);
  v = -cos(x) * sin(y) * cos(z);
  w = 0.0;
}

// Set BC type on each boundary face
DEVICE_PREFIX void ps3d_set_boundary_type(const DG_FP x0, const DG_FP y0, const DG_FP z0,
                                          const DG_FP x1, const DG_FP y1, const DG_FP z1,
                                          const DG_FP x2, const DG_FP y2, const DG_FP z2,
                                          int &bc_type) {
  if(fp_equal(x0, x1) && fp_equal(x0, x2) && x0 < -LW_INLET_LENGTH + 0.1) {
    bc_type = BC_TYPE_INFLOW;
  } else if(x0 < 1e-8 && x1 < 1e-8 && x2 < 1e-8) {
    bc_type = BC_TYPE_NO_SLIP;
  } else if(fp_equal(x0, x1) && fp_equal(x0, x2) && x0 > LW_LENGTH_SHORT - 0.1) {
    bc_type = BC_TYPE_NATURAL_OUTFLOW;
  } else if(y0 * y0 + z0 * z0 < LW_RADIUS * LW_RADIUS - 5e-1 &&
            y1 * y1 + z1 * z1 < LW_RADIUS * LW_RADIUS - 5e-1 &&
            y2 * y2 + z2 * z2 < LW_RADIUS * LW_RADIUS - 5e-1) {
    bc_type = BC_TYPE_NO_SLIP;
  } else {
    bc_type = BC_TYPE_SLIP;
  }
}

// Custom BC pressure and viscosity linear solves BC conditions
DEVICE_PREFIX void ps3d_custom_bc_get_pr_type(const int bc_type, int &pr_bc) {
  if(bc_type == BC_TYPE_INFLOW) {
    pr_bc = BC_NEUMANN;
  }
}

DEVICE_PREFIX void ps3d_custom_bc_get_vis_type(const int bc_type, int &vis_bc) {
  if(bc_type == BC_TYPE_INFLOW) {
    vis_bc = BC_DIRICHLET;
  }
}

// Custom BC velocities on boundary
DEVICE_PREFIX void ps3d_custom_bc_get_vel(const int bc_type, const DG_FP time,
                                          const DG_FP x, const DG_FP y, const DG_FP z,
                                          const DG_FP nx, const DG_FP ny, const DG_FP nz,
                                          const DG_FP mU, const DG_FP mV, const DG_FP mW,
                                          DG_FP &u, DG_FP &v, DG_FP &w) {
  if(bc_type == BC_TYPE_INFLOW) {
    DG_FP tmp = y * y + z * z - LW_INLET_RADIUS * LW_INLET_RADIUS;
    tmp = fabs(tmp);
    tmp = fmin(1.0, tmp / 0.1);
    u = tmp * 1.0;
    v = 0.0;
    w = 0.0;
  }
}

// Custom BC pressure Neumann boundary condition (return 0.0 if not Neumann) [SINGLE-PHASE]
DEVICE_PREFIX DG_FP ps3d_custom_bc_get_pr_neumann(const int bc_type, const DG_FP time,
                          const DG_FP x, const DG_FP y, const DG_FP z, const DG_FP nx,
                          const DG_FP ny, const DG_FP nz, const DG_FP N0, const DG_FP N1,
                          const DG_FP N2, const DG_FP curl20, const DG_FP curl21,
                          const DG_FP curl22, const DG_FP reynolds) {
  if(bc_type == BC_TYPE_INFLOW) {
    DG_FP res0 = -N0 - curl20 / reynolds;
    DG_FP res1 = -N1 - curl21 / reynolds;
    DG_FP res2 = -N2 - curl22 / reynolds;
    return nx * res0 + ny * res1 + nz * res2;
  }

  return 0.0;
}

// Custom BC pressure Dirchlet boundary condition (return 0.0 if not Dirchlet)
DEVICE_PREFIX DG_FP ps3d_custom_bc_get_pr_dirichlet(const int bc_type, const DG_FP time,
                                                    const DG_FP x, const DG_FP y, const DG_FP z) {
  return 0.0;
}

// Custom BC viscosity Neumann boundary conditions (return 0.0 if not Neumann)
DEVICE_PREFIX void ps3d_custom_bc_get_vis_neumann(const int bc_type, const DG_FP time,
                                          const DG_FP x, const DG_FP y, const DG_FP z,
                                          const DG_FP nx, const DG_FP ny, const DG_FP nz,
                                          const DG_FP mU, const DG_FP mV, const DG_FP mW,
                                          DG_FP &u, DG_FP &v, DG_FP &w) {
  u = 0.0;
  v = 0.0;
  w = 0.0;
}

// Set the initial interface between phases for multiphase simulations
DEVICE_PREFIX void ps3d_set_surface(const DG_FP x, const DG_FP y, const DG_FP z,
                                    DG_FP &s) {
  // s[i] = sqrt((x[i] - 3.0) * (x[i] - 3.0) + (y[i] - 3.0) * (y[i] - 3.0) + (z[i] - 3.0) * (z[i] - 3.0)) - 1.5;
  // s[i] = x[i] - 1.0;
  // s[i] = sqrt((x[i] - 1.0) * (x[i] - 1.0) + (y[i] - 0.5) * (y[i] - 0.5) + (z[i] - 0.5) * (z[i] - 0.5)) - 0.25;
  // s[i] = sqrt((x[i] - 0.1) * (x[i] - 0.1) + y[i] * y[i] + z[i] * z[i]) - 0.05;
  // if(x[i] >= -1e-5)
  //   s[i] = fmax(fmax(fabs(x[i] + 1e-3), fabs(y[i])), fabs(z[i]));
  // else
  //   s[i] = x[i] + 1e-3;

  s = sqrt((x + 3.0) * (x + 3.0) + y * y + z * z) - 2.99;
  s = fmax(fmin(1.0, s), -1.0);
}

// Set level set value on custom BCs (return sM otherwise)
DEVICE_PREFIX DG_FP ps3d_custom_bc_get_ls(const int bc_type, const DG_FP x,
                                const DG_FP y, const DG_FP z, const DG_FP sM) {
  if(bc_type == BC_TYPE_INFLOW) {
    return -1.0;
  }

  return sM;
}

// Custom BC pressure Neumann boundary condition (return 0.0 if not Neumann) [MULTI-PHASE]
DEVICE_PREFIX DG_FP ps3d_custom_bc_get_pr_neumann_multiphase(const int bc_type, const DG_FP time,
                          const DG_FP x, const DG_FP y, const DG_FP z, const DG_FP nx,
                          const DG_FP ny, const DG_FP nz, const DG_FP N0, const DG_FP N1,
                          const DG_FP N2, const DG_FP curl20, const DG_FP curl21,
                          const DG_FP curl22, const DG_FP reynolds, const DG_FP rho) {
  if(bc_type == BC_TYPE_INFLOW) {
    DG_FP res0 = -N0 - curl20 / (reynolds * rho);
    DG_FP res1 = -N1 - curl21 / (reynolds * rho);
    DG_FP res2 = -N2 - curl22 / (reynolds * rho);
    return nx * res0 + ny * res1 + nz * res2;
  }

  return 0.0;
}

#endif
#endif
