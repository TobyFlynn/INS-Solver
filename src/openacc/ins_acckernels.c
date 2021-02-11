//
// auto-generated by op2.py
//

// global constants
extern double gam;
extern double mu;
extern double bc_mach;
extern double bc_alpha;
extern double bc_p;
extern double bc_r;
extern double bc_u;
extern double bc_v;
extern double bc_e;
extern double ones[15];
extern int FMASK[15];

// header
#include "op_lib_c.h"

void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){}
// user kernel files
#include "init_grid_acckernel.c"
#include "set_ic_acckernel.c"
