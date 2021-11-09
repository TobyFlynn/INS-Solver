//
// auto-generated by op2.py
//

//user function
#include "../kernels/zero_npf1.h"

// host stub function
void op_par_loop_zero_npf1(char const *name, op_set set,
  op_arg arg0){

  int nargs = 1;
  op_arg args[1];

  args[0] = arg0;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(34);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  zero_npf1");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      zero_npf1(
        &((double*)arg0.data)[12*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[34].name      = name;
  OP_kernels[34].count    += 1;
  OP_kernels[34].time     += wall_t2 - wall_t1;
  OP_kernels[34].transfer += (float)set->size * arg0.size * 2.0f;
}
