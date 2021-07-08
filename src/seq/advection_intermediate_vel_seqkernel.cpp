//
// auto-generated by op2.py
//

//user function
#include "../kernels/advection_intermediate_vel.h"

// host stub function
void op_par_loop_advection_intermediate_vel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15){

  int nargs = 16;
  op_arg args[16];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(32);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  advection_intermediate_vel");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      advection_intermediate_vel(
        (double*)arg0.data,
        (double*)arg1.data,
        (double*)arg2.data,
        (double*)arg3.data,
        (double*)arg4.data,
        (double*)arg5.data,
        &((double*)arg6.data)[15*n],
        &((double*)arg7.data)[15*n],
        &((double*)arg8.data)[15*n],
        &((double*)arg9.data)[15*n],
        &((double*)arg10.data)[15*n],
        &((double*)arg11.data)[15*n],
        &((double*)arg12.data)[15*n],
        &((double*)arg13.data)[15*n],
        &((double*)arg14.data)[15*n],
        &((double*)arg15.data)[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[32].name      = name;
  OP_kernels[32].count    += 1;
  OP_kernels[32].time     += wall_t2 - wall_t1;
  OP_kernels[32].transfer += (float)set->size * arg6.size;
  OP_kernels[32].transfer += (float)set->size * arg7.size;
  OP_kernels[32].transfer += (float)set->size * arg8.size;
  OP_kernels[32].transfer += (float)set->size * arg9.size;
  OP_kernels[32].transfer += (float)set->size * arg10.size;
  OP_kernels[32].transfer += (float)set->size * arg11.size;
  OP_kernels[32].transfer += (float)set->size * arg12.size;
  OP_kernels[32].transfer += (float)set->size * arg13.size;
  OP_kernels[32].transfer += (float)set->size * arg14.size * 2.0f;
  OP_kernels[32].transfer += (float)set->size * arg15.size * 2.0f;
}
