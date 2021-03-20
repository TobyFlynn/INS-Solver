//
// auto-generated by op2.py
//

//user function
#include "../kernels/gauss_tau_bc.h"

// host stub function
void op_par_loop_gauss_tau_bc(char const *name, op_set set,
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
  op_timing_realloc(22);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_tau_bc\n");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map1idx;
      map1idx = arg1.map_data[n * arg1.map->dim + 0];


      gauss_tau_bc(
        &((int*)arg0.data)[1 * n],
        &((double*)arg1.data)[15 * map1idx],
        &((double*)arg2.data)[3 * map1idx]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[22].name      = name;
  OP_kernels[22].count    += 1;
  OP_kernels[22].time     += wall_t2 - wall_t1;
  OP_kernels[22].transfer += (float)set->size * arg1.size;
  OP_kernels[22].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[22].transfer += (float)set->size * arg0.size;
  OP_kernels[22].transfer += (float)set->size * arg1.map->dim * 4.0f;
}
