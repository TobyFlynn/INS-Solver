#ifndef __INS_PROBLEM_SPECIFIC_3D_H
#define __INS_PROBLEM_SPECIFIC_3D_H

#if defined(INS_CUDA)
#define DEVICE_PREFIX __device__
#elif defined(INS_HIP)
#define DEVICE_PREFIX __device__
#else
#define DEVICE_PREFIX inline
#endif

#define LW_INFLOW_BC 0
#define LW_OUTFLOW_BC 1
#define LW_SLIP_WALL_BC 2
#define LW_NO_SLIP_WALL_BC 3

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

#endif
