#ifndef __INS_3D_SCRUBBER_H
#define __INS_3D_SCRUBBER_H

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
#define LS_CAP 0.25
// Problem specifc definitions
#define DOMAIN_HEIGHT 4.203912

/************************************************************************
 * You can edit the body of the functions below but not their signature *
 ************************************************************************/

// Set the initial conditions of the problem
DEVICE_PREFIX void ps3d_set_ic(const DG_FP x, const DG_FP y, const DG_FP z,
                               DG_FP &u, DG_FP &v, DG_FP &w) {
  // Liquid Whistle Initial Conditions
  u = 0.0;
  v = 0.0;
  w = 0.0;
}

// Set BC type on each boundary face
DEVICE_PREFIX void ps3d_set_boundary_type(const DG_FP x0, const DG_FP y0, const DG_FP z0,
                                          const DG_FP x1, const DG_FP y1, const DG_FP z1,
                                          const DG_FP x2, const DG_FP y2, const DG_FP z2,
                                          int &bc_type) {
  bc_type = BC_TYPE_SLIP;
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
  const DG_FP MIN_X = -1.19;
  const DG_FP MAX_X = 1.18;
  const DG_FP IC_DEPTH = 0.15;
  // Plate 0
  const DG_FP MIN_Y_PLATE_0 = 1.62;
  const DG_FP MAX_Y_PLATE_0 = MIN_Y_PLATE_0 + IC_DEPTH;
  // Plate 1
  const DG_FP MIN_Y_PLATE_1 = 0.55;
  const DG_FP MAX_Y_PLATE_1 = MIN_Y_PLATE_1 + IC_DEPTH;
  // Plate 2
  const DG_FP MIN_Y_PLATE_2 = -0.47;
  const DG_FP MAX_Y_PLATE_2 = MIN_Y_PLATE_2 + IC_DEPTH;
  // Plate 3
  const DG_FP MIN_Y_PLATE_3 = -1.54;
  const DG_FP MAX_Y_PLATE_3 = MIN_Y_PLATE_3 + IC_DEPTH;

  if(x > MIN_X && x < MAX_X) {
    const DG_FP min_dist_to_x_boundary = fmin(x - MIN_X, MAX_X - x);
    if(y > MIN_Y_PLATE_0 && y < MAX_Y_PLATE_0) {
      s = -fmin(MAX_Y_PLATE_0 - y, min_dist_to_x_boundary);
    } else if(y > MIN_Y_PLATE_1 && y < MAX_Y_PLATE_1) {
      s = -fmin(MAX_Y_PLATE_1 - y, min_dist_to_x_boundary);
    } else if(y > MIN_Y_PLATE_2 && y < MAX_Y_PLATE_2) {
      s = -fmin(MAX_Y_PLATE_2 - y, min_dist_to_x_boundary);
    } else if(y > MIN_Y_PLATE_3 && y < MAX_Y_PLATE_3) {
      s = -fmin(MAX_Y_PLATE_3 - y, min_dist_to_x_boundary);
    } else {
      if(y > MAX_Y_PLATE_0) {
        s = y - MAX_Y_PLATE_0;
      } else if(y > MAX_Y_PLATE_1) {
        s = y - MAX_Y_PLATE_1;
      } else if(y > MAX_Y_PLATE_2) {
        s = y - MAX_Y_PLATE_2;
      } else if(y > MAX_Y_PLATE_3) {
        s = y - MAX_Y_PLATE_3;
      } else {
        s = (y + DOMAIN_HEIGHT) - MAX_Y_PLATE_0;
      }
    }
  } else {
    if(y > MIN_Y_PLATE_0 && y < MAX_Y_PLATE_0 || y > MIN_Y_PLATE_1 && y < MAX_Y_PLATE_1
       || y > MIN_Y_PLATE_2 && y < MAX_Y_PLATE_2 || y > MIN_Y_PLATE_3 && y < MAX_Y_PLATE_3) {
      s = x < MIN_X ? MIN_X - x : x - MAX_X;
    } else {
      if(y > MAX_Y_PLATE_0) {
        DG_FP dist_pt_l = (x - MIN_X) * (x - MIN_X) + (y - MAX_Y_PLATE_0) * (y - MAX_Y_PLATE_0);
        DG_FP dist_pt_r = (x - MAX_X) * (x - MAX_X) + (y - MAX_Y_PLATE_0) * (y - MAX_Y_PLATE_0);
        s = sqrt(fmin(dist_pt_l, dist_pt_r));
      } else if(y > MAX_Y_PLATE_1) {
        DG_FP dist_pt_l = (x - MIN_X) * (x - MIN_X) + (y - MAX_Y_PLATE_1) * (y - MAX_Y_PLATE_1);
        DG_FP dist_pt_r = (x - MAX_X) * (x - MAX_X) + (y - MAX_Y_PLATE_1) * (y - MAX_Y_PLATE_1);
        s = sqrt(fmin(dist_pt_l, dist_pt_r));
      } else if(y > MAX_Y_PLATE_2) {
        DG_FP dist_pt_l = (x - MIN_X) * (x - MIN_X) + (y - MAX_Y_PLATE_2) * (y - MAX_Y_PLATE_2);
        DG_FP dist_pt_r = (x - MAX_X) * (x - MAX_X) + (y - MAX_Y_PLATE_2) * (y - MAX_Y_PLATE_2);
        s = sqrt(fmin(dist_pt_l, dist_pt_r));
      } else if(y > MAX_Y_PLATE_3) {
        DG_FP dist_pt_l = (x - MIN_X) * (x - MIN_X) + (y - MAX_Y_PLATE_3) * (y - MAX_Y_PLATE_3);
        DG_FP dist_pt_r = (x - MAX_X) * (x - MAX_X) + (y - MAX_Y_PLATE_3) * (y - MAX_Y_PLATE_3);
        s = sqrt(fmin(dist_pt_l, dist_pt_r));
      } else {
        DG_FP dist_pt_l = (x - MIN_X) * (x - MIN_X) + ((y + DOMAIN_HEIGHT) - MAX_Y_PLATE_0) * ((y + DOMAIN_HEIGHT) - MAX_Y_PLATE_0);
        DG_FP dist_pt_r = (x - MAX_X) * (x - MAX_X) + ((y + DOMAIN_HEIGHT) - MAX_Y_PLATE_0) * ((y + DOMAIN_HEIGHT) - MAX_Y_PLATE_0);
        s = sqrt(fmin(dist_pt_l, dist_pt_r));
      }
    }
  }
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
  u = 0.0;
  v = 0.0;
  w = 0.0;
}

#endif
