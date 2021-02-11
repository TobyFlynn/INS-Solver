//
// auto-generated by op2.py
//

#ifdef _OPENMP
  #include <omp.h>
#endif

// global constants
extern double gam;
extern double mu;
extern double bc_mach;
extern double bc_alpha;
extern double bc_p;
extern double bc_u;
extern double bc_v;
extern double ones[15];
extern int FMASK[15];

// header
#include "op_lib_cpp.h"

// user kernel files
#include "init_grid_kernel.cpp"
#include "set_ic_kernel.cpp"
#include "advection_flux_kernel.cpp"
