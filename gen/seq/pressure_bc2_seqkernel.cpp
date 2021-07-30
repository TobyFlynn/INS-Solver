//
// auto-generated by op2.py
//

//user function
#include "../kernels/pressure_bc2.h"

// host stub function
void op_par_loop_pressure_bc2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  int nargs = 8;
  op_arg args[8];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(35);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_bc2\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map4idx;
      map4idx = arg4.map_data[n * arg4.map->dim + 0];


      pressure_bc2(
        &((int*)arg0.data)[1 * n],
        &((int*)arg1.data)[1 * n],
        (double*)arg2.data,
        (int*)arg3.data,
        &((double*)arg4.data)[12 * map4idx],
        &((double*)arg5.data)[12 * map4idx],
        &((double*)arg6.data)[12 * map4idx],
        &((double*)arg7.data)[12 * map4idx]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[35].name      = name;
  OP_kernels[35].count    += 1;
  OP_kernels[35].time     += wall_t2 - wall_t1;
  OP_kernels[35].transfer += (float)set->size * arg4.size;
  OP_kernels[35].transfer += (float)set->size * arg5.size;
  OP_kernels[35].transfer += (float)set->size * arg6.size;
  OP_kernels[35].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[35].transfer += (float)set->size * arg0.size;
  OP_kernels[35].transfer += (float)set->size * arg1.size;
  OP_kernels[35].transfer += (float)set->size * arg2.size;
  OP_kernels[35].transfer += (float)set->size * arg3.size;
  OP_kernels[35].transfer += (float)set->size * arg4.map->dim * 4.0f;
}
