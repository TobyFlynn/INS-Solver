#ifndef __INS_PROBLEM_SPECIFIC_2D_H
#define __INS_PROBLEM_SPECIFIC_2D_H

#if defined(INS_CUDA)
#define DEVICE_PREFIX __device__
#elif defined(INS_HIP)
#define DEVICE_PREFIX __device__
#else
#define DEVICE_PREFIX inline
#endif

/************************************************************************
 * You can edit the body of the functions below but not their signature *
 ************************************************************************/

// Set the initial conditions of the problem
DEVICE_PREFIX void ps2d_set_ic(const DG_FP x, const DG_FP y, DG_FP &u, DG_FP &v) {
  // Euler vortex initial conditions
  // const DG_FP PI = 3.141592653589793238463;
  // const DG_FP R = 1.5;
  // const DG_FP S = 13.5;
  // DG_FP f = (1.0 - x * x - y * y) / (2.0 * R * R);
  // u = (S * y * exp(f)) / (2.0 * PI * R);
  // v = (1.0 - ((S * x * exp(f)) / (2.0 * PI * R)));

  u = 0.0;
  v = 0.0;
}

// Set the initial interface between phases for multiphase simulations
DEVICE_PREFIX void ps2d_set_surface(const DG_FP x, const DG_FP y, DG_FP &s) {
  s = sqrt(x * x + y * y) - 7.5;
}

#endif
