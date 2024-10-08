#ifndef __INS_3D_ZALESAK_H
#define __INS_3D_ZALESAK_H

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
#define BC_SLIP 2

// Hardcoded BC types, do not edit
#define BC_TYPE_NO_SLIP 0
#define BC_TYPE_SLIP 1
#define BC_TYPE_NATURAL_OUTFLOW 2

// Required definitions
#define LS_CAP 0.1
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
  const DG_FP xO = 0.5;
  const DG_FP yO = 0.75;
  const DG_FP r = 0.15;
  const DG_FP w = 0.05;
  const DG_FP h = 0.25;

  const DG_FP xL = xO - 0.5*w;
  const DG_FP xR = xO + 0.5*w;

  const DG_FP yT = yO - r + h;
  const DG_FP yB = yO - sqrt(r*r - 0.5*w*w);

  DG_FP f1 = sqrt((x - xO)*(x - xO) + (y - yO)*(y - yO)) - r;
  DG_FP f2 = x - xL;
  DG_FP f3 = -(x - xR);
  DG_FP f4 = -(y - yT);
  DG_FP f5 = sqrt((x - xL)*(x - xL) + (y - yB)*(y - yB));
  DG_FP f6 = sqrt((x - xR)*(x - xR) + (y - yB)*(y - yB));
  DG_FP f7 = -sqrt((x - xL)*(x - xL) + (y - yT)*(y - yT));
  DG_FP f8 = -sqrt((x - xR)*(x - xR) + (y - yT)*(y - yT));

  if (y >= yT)
  {
      if (x <= xL)
          s = fabs(f1) < fabs(f7) ? f1 : f7;
      else if (x > xR)
          s = fabs(f1) < fabs(f8) ? f1 : f8;
      else
          s = fabs(f1) < fabs(f4) ? f1 : f4;
  }
  else if ((y < yT) && (y >= yB))
  {
      if (x <= xL)
          s = fabs(f1) < fabs(f2) ? f1 : f2;
      else if (x > xR)
          s = fabs(f1) < fabs(f3) ? f1 : f3;
      else
      {
          if ((fabs(f4) < fabs(f2)) && (fabs(f4) < fabs(f3)))
              s = f4;
          else
              s = fabs(f2) < fabs(f3) ? f2 : f3;
      }
  }
  else if (y < yB)
  {
      if ((x <= xL) || (x > xR))
          s = f1;
      else
          s = fabs(f5) < fabs(f6) ? f5 : f6;
  }
  else
      s = 0.0;

  s = fmax(fmin(LS_CAP, s), -LS_CAP);
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

// Set the velocity for the level-set-only solver based on current time and coordinates
DEVICE_PREFIX void ps3d_set_ls_vel(const DG_FP time, const DG_FP x, const DG_FP y, const DG_FP z, 
                                   DG_FP &u, DG_FP &v, DG_FP &w) {
  const DG_FP PI = 3.141592653589793238463;
  u = 2.0 * PI * (0.5 - y);
  v = 2.0 * PI * (x - 0.5);
  w = 0.0;
}

#endif
