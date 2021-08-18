//
// auto-generated by op2.py
//

//user function
#include "../kernels/ls_normalise.h"

// host stub function
void op_par_loop_ls_normalise(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(70);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_normalise");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      ls_normalise(
        &((double*)arg0.data)[10*n],
        &((double*)arg1.data)[10*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[70].name      = name;
  OP_kernels[70].count    += 1;
  OP_kernels[70].time     += wall_t2 - wall_t1;
  OP_kernels[70].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[70].transfer += (float)set->size * arg1.size * 2.0f;
}