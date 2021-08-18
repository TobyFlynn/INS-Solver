//
// auto-generated by op2.py
//

//user function
#include "../kernels/ls_advec_bedges.h"

// host stub function
void op_par_loop_ls_advec_bedges(char const *name, op_set set,
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
  op_timing_realloc(54);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_advec_bedges\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map2idx;
      map2idx = arg2.map_data[n * arg2.map->dim + 0];


      ls_advec_bedges(
        &((int*)arg0.data)[1 * n],
        &((int*)arg1.data)[1 * n],
        &((double*)arg2.data)[10 * map2idx],
        &((double*)arg3.data)[10 * map2idx],
        &((double*)arg4.data)[10 * map2idx],
        &((double*)arg5.data)[12 * map2idx]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[54].name      = name;
  OP_kernels[54].count    += 1;
  OP_kernels[54].time     += wall_t2 - wall_t1;
  OP_kernels[54].transfer += (float)set->size * arg2.size;
  OP_kernels[54].transfer += (float)set->size * arg3.size;
  OP_kernels[54].transfer += (float)set->size * arg4.size;
  OP_kernels[54].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[54].transfer += (float)set->size * arg0.size;
  OP_kernels[54].transfer += (float)set->size * arg1.size;
  OP_kernels[54].transfer += (float)set->size * arg2.map->dim * 4.0f;
}