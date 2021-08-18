//
// auto-generated by op2.py
//

//user function
#include "../kernels/poisson_edges.h"

// host stub function
void op_par_loop_poisson_edges(char const *name, op_set set,
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
  op_timing_realloc(20);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_edges\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map0idx;
      int map3idx;
      map0idx = arg0.map_data[n * arg0.map->dim + 0];
      map3idx = arg0.map_data[n * arg0.map->dim + 1];


      poisson_edges(
        &((double*)arg0.data)[10 * map0idx],
        &((double*)arg1.data)[100 * n],
        &((double*)arg2.data)[10 * map0idx],
        &((double*)arg0.data)[10 * map3idx],
        &((double*)arg4.data)[100 * n],
        &((double*)arg2.data)[10 * map3idx]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[20].name      = name;
  OP_kernels[20].count    += 1;
  OP_kernels[20].time     += wall_t2 - wall_t1;
  OP_kernels[20].transfer += (float)set->size * arg0.size;
  OP_kernels[20].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg1.size;
  OP_kernels[20].transfer += (float)set->size * arg4.size;
  OP_kernels[20].transfer += (float)set->size * arg0.map->dim * 4.0f;
}
