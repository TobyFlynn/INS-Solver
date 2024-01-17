#ifndef __INS_3D_LS_TEST_H
#define __INS_3D_LS_TEST_H

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

// Required definitions
#define LS_CAP 100.0
// Problem specifc definitions

/************************************************************************
 * You can edit the body of the functions below but not their signature *
 ************************************************************************/

// Set the initial conditions of the problem
DEVICE_PREFIX void ps3d_set_ic(const DG_FP x, const DG_FP y, const DG_FP z,
                               DG_FP &u, DG_FP &v, DG_FP &w) {
  u = 0.0;
  v = 0.0;
  w = 0.0;
}

// Set BC type on each boundary face
DEVICE_PREFIX void ps3d_set_boundary_type(const DG_FP x0, const DG_FP y0, const DG_FP z0,
                                          const DG_FP x1, const DG_FP y1, const DG_FP z1,
                                          const DG_FP x2, const DG_FP y2, const DG_FP z2,
                                          int &bc_type) {
  bc_type = BC_TYPE_NATURAL_OUTFLOW;
}

// Custom BC pressure and viscosity linear solves BC conditions
DEVICE_PREFIX void ps3d_custom_bc_get_pr_type(const int bc_type, int &pr_bc) {

}

DEVICE_PREFIX void ps3d_custom_bc_get_vis_type(const int bc_type, int &vis_bc) {

}

// Custom BC velocities on boundary
DEVICE_PREFIX void ps3d_custom_bc_get_vel(const int bc_type, const DG_FP time,
                                          const DG_FP x, const DG_FP y, const DG_FP z,
                                          const DG_FP nx, const DG_FP ny, const DG_FP nz,
                                          const DG_FP mU, const DG_FP mV, const DG_FP mW,
                                          DG_FP &u, DG_FP &v, DG_FP &w) {

}

// Custom BC pressure Neumann boundary condition (return 0.0 if not Neumann) [SINGLE-PHASE]
DEVICE_PREFIX DG_FP ps3d_custom_bc_get_pr_neumann(const int bc_type, const DG_FP time,
                          const DG_FP x, const DG_FP y, const DG_FP z, const DG_FP nx,
                          const DG_FP ny, const DG_FP nz, const DG_FP N0, const DG_FP N1,
                          const DG_FP N2, const DG_FP curl20, const DG_FP curl21,
                          const DG_FP curl22, const DG_FP reynolds) {
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
  s = sqrt((x - 0.5) * (x - 0.5) + (y - 0.75) * (y - 0.75)) - 0.15;
}

// Set level set value on custom BCs (return sM otherwise)
DEVICE_PREFIX DG_FP ps3d_custom_bc_get_ls(const int bc_type, const DG_FP x,
                                const DG_FP y, const DG_FP z, const DG_FP sM) {
  return sM;
}

// Custom BC pressure Neumann boundary condition (return 0.0 if not Neumann) [MULTI-PHASE]
DEVICE_PREFIX DG_FP ps3d_custom_bc_get_pr_neumann_multiphase(const int bc_type, const DG_FP time,
                          const DG_FP x, const DG_FP y, const DG_FP z, const DG_FP nx,
                          const DG_FP ny, const DG_FP nz, const DG_FP N0, const DG_FP N1,
                          const DG_FP N2, const DG_FP curl20, const DG_FP curl21,
                          const DG_FP curl22, const DG_FP reynolds, const DG_FP rho) {
  return 0.0;
}

#endif
