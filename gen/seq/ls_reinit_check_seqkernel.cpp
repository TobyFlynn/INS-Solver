//
// auto-generated by op2.py
//

//user function
#include "../kernels/ls_reinit_check.h"

// host stub function
void op_par_loop_ls_reinit_check(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(68);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_reinit_check");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      ls_reinit_check(
        (double*)arg0.data,
        &((double*)arg1.data)[10*n],
        &((double*)arg2.data)[10*n],
        &((double*)arg3.data)[10*n],
        (double*)arg4.data,
        (int*)arg5.data);
    }
  }

  // combine reduction data
  op_mpi_reduce_double(&arg4,(double*)arg4.data);
  op_mpi_reduce_int(&arg5,(int*)arg5.data);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[68].name      = name;
  OP_kernels[68].count    += 1;
  OP_kernels[68].time     += wall_t2 - wall_t1;
  OP_kernels[68].transfer += (float)set->size * arg1.size;
  OP_kernels[68].transfer += (float)set->size * arg2.size;
  OP_kernels[68].transfer += (float)set->size * arg3.size;
}
