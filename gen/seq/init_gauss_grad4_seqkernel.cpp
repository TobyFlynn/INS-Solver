//
// auto-generated by op2.py
//

//user function
#include "../kernels/init_gauss_grad4.h"

// host stub function
void op_par_loop_init_gauss_grad4(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  int nargs = 10;
  op_arg args[10];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(7);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: init_gauss_grad4\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map1idx;
      map1idx = arg1.map_data[n * arg1.map->dim + 0];


      init_gauss_grad4(
        &((int*)arg0.data)[1 * n],
        &((double*)arg1.data)[18 * map1idx],
        &((double*)arg2.data)[18 * map1idx],
        &((double*)arg3.data)[60 * map1idx],
        &((double*)arg4.data)[60 * map1idx],
        &((double*)arg5.data)[60 * map1idx],
        &((double*)arg6.data)[60 * map1idx],
        &((double*)arg7.data)[60 * map1idx],
        &((double*)arg8.data)[60 * map1idx],
        &((double*)arg9.data)[60 * n]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[7].name      = name;
  OP_kernels[7].count    += 1;
  OP_kernels[7].time     += wall_t2 - wall_t1;
  OP_kernels[7].transfer += (float)set->size * arg1.size;
  OP_kernels[7].transfer += (float)set->size * arg2.size;
  OP_kernels[7].transfer += (float)set->size * arg3.size;
  OP_kernels[7].transfer += (float)set->size * arg4.size;
  OP_kernels[7].transfer += (float)set->size * arg5.size;
  OP_kernels[7].transfer += (float)set->size * arg6.size;
  OP_kernels[7].transfer += (float)set->size * arg7.size;
  OP_kernels[7].transfer += (float)set->size * arg8.size;
  OP_kernels[7].transfer += (float)set->size * arg0.size;
  OP_kernels[7].transfer += (float)set->size * arg9.size;
  OP_kernels[7].transfer += (float)set->size * arg1.map->dim * 4.0f;
}
