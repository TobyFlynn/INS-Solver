//
// auto-generated by op2.py
//

//user function
#include "../kernels/init_cubature_grad.h"

// host stub function
void op_par_loop_init_cubature_grad(char const *name, op_set set,
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
  op_timing_realloc(1);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_cubature_grad");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      init_cubature_grad(
        &((double*)arg0.data)[16*n],
        &((double*)arg1.data)[16*n],
        &((double*)arg2.data)[16*n],
        &((double*)arg3.data)[16*n],
        &((double*)arg4.data)[96*n],
        &((double*)arg5.data)[96*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[1].name      = name;
  OP_kernels[1].count    += 1;
  OP_kernels[1].time     += wall_t2 - wall_t1;
  OP_kernels[1].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[1].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[1].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[1].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[1].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[1].transfer += (float)set->size * arg5.size * 2.0f;
}
