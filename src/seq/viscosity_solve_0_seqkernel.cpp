//
// auto-generated by op2.py
//

//user function
#include "../kernels/viscosity_solve_0.h"

// host stub function
void op_par_loop_viscosity_solve_0(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  int nargs = 9;
  op_arg args[9];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(31);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  viscosity_solve_0");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      viscosity_solve_0(
        &((double*)arg0.data)[46*n],
        &((double*)arg1.data)[690*n],
        &((double*)arg2.data)[690*n],
        &((double*)arg3.data)[46*n],
        &((double*)arg4.data)[46*n],
        (double*)arg5.data,
        &((double*)arg6.data)[225*n],
        &((double*)arg7.data)[15*n],
        &((double*)arg8.data)[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[31].name      = name;
  OP_kernels[31].count    += 1;
  OP_kernels[31].time     += wall_t2 - wall_t1;
  OP_kernels[31].transfer += (float)set->size * arg0.size;
  OP_kernels[31].transfer += (float)set->size * arg1.size;
  OP_kernels[31].transfer += (float)set->size * arg2.size;
  OP_kernels[31].transfer += (float)set->size * arg3.size;
  OP_kernels[31].transfer += (float)set->size * arg4.size;
  OP_kernels[31].transfer += (float)set->size * arg6.size;
  OP_kernels[31].transfer += (float)set->size * arg7.size;
  OP_kernels[31].transfer += (float)set->size * arg8.size * 2.0f;
}
