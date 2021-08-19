//
// auto-generated by op2.py
//

//user function
#include "../kernels/glb_ind_kernel.h"

// host stub function
void op_par_loop_glb_ind_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  arg0.idx = 0;
  args[0] = arg0;
  for ( int v=1; v<2; v++ ){
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 1, "int", OP_READ);
  }

  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(9);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: glb_ind_kernel\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map0idx;
      int map1idx;
      map0idx = arg0.map_data[n * arg0.map->dim + 0];
      map1idx = arg0.map_data[n * arg0.map->dim + 1];

      const int* arg0_vec[] = {
         &((int*)arg0.data)[1 * map0idx],
         &((int*)arg0.data)[1 * map1idx]};

      glb_ind_kernel(
        arg0_vec,
        &((int*)arg2.data)[1 * n],
        &((int*)arg3.data)[1 * n]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[9].name      = name;
  OP_kernels[9].count    += 1;
  OP_kernels[9].time     += wall_t2 - wall_t1;
  OP_kernels[9].transfer += (float)set->size * arg0.size;
  OP_kernels[9].transfer += (float)set->size * arg2.size;
  OP_kernels[9].transfer += (float)set->size * arg3.size;
  OP_kernels[9].transfer += (float)set->size * arg0.map->dim * 4.0f;
}
