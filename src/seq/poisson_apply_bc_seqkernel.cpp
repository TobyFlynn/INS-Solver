//
// auto-generated by op2.py
//

//user function
#include "../kernels/poisson_apply_bc.h"

// host stub function
void op_par_loop_poisson_apply_bc(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(26);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_apply_bc\n");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map2idx;
      map2idx = arg2.map_data[n * arg2.map->dim + 0];


      poisson_apply_bc(
        &((int*)arg0.data)[1 * n],
        &((double*)arg1.data)[105 * n],
        &((double*)arg2.data)[21 * map2idx],
        &((double*)arg3.data)[15 * map2idx]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[26].name      = name;
  OP_kernels[26].count    += 1;
  OP_kernels[26].time     += wall_t2 - wall_t1;
  OP_kernels[26].transfer += (float)set->size * arg2.size;
  OP_kernels[26].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[26].transfer += (float)set->size * arg0.size;
  OP_kernels[26].transfer += (float)set->size * arg1.size;
  OP_kernels[26].transfer += (float)set->size * arg2.map->dim * 4.0f;
}
