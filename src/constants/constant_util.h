#ifndef __CONSTANTS_UTIL_H
#define __CONSTANTS_UTIL_H

// Physics constants
double gam, mu, nu0, nu1, rho0, rho1, bc_mach, bc_alpha, bc_p, bc_r, bc_u, bc_v, bc_e, dt;
double ic_u, ic_v;
double vortex_x0 = 5.0;
double vortex_y0 = 0.0;
double vortex_beta = 5.0;

double lift_drag_vec[5] = {0.10000000000000017, 0.5444444444444444, 0.711111111111111, 0.54444444444444462, 0.10000000000000016};

#endif
