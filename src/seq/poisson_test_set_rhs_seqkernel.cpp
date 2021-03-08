//
// auto-generated by op2.py
//

//user function
#include "../kernels/poisson_test_set_rhs.h"

// host stub function
void op_par_loop_poisson_test_set_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(32);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_test_set_rhs");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      poisson_test_set_rhs(
        &((double*)arg0.data)[15*n],
        &((double*)arg1.data)[15*n],
        &((double*)arg2.data)[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[32].name      = name;
  OP_kernels[32].count    += 1;
  OP_kernels[32].time     += wall_t2 - wall_t1;
  OP_kernels[32].transfer += (float)set->size * arg0.size;
  OP_kernels[32].transfer += (float)set->size * arg1.size;
  OP_kernels[32].transfer += (float)set->size * arg2.size * 2.0f;
}
