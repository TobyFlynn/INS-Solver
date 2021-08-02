//
// auto-generated by op2.py
//

//user function
#include "../kernels/gauss_tau.h"

// host stub function
void op_par_loop_gauss_tau(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  arg1.idx = 0;
  args[1] = arg1;
  for ( int v=1; v<2; v++ ){
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 12, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 3, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(4);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_tau\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map1idx;
      int map2idx;
      map1idx = arg1.map_data[n * arg1.map->dim + 0];
      map2idx = arg1.map_data[n * arg1.map->dim + 1];

      const double* arg1_vec[] = {
         &((double*)arg1.data)[12 * map1idx],
         &((double*)arg1.data)[12 * map2idx]};
      double* arg3_vec[] = {
         &((double*)arg3.data)[3 * map1idx],
         &((double*)arg3.data)[3 * map2idx]};

      gauss_tau(
        &((int*)arg0.data)[2 * n],
        arg1_vec,
        arg3_vec);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[4].name      = name;
  OP_kernels[4].count    += 1;
  OP_kernels[4].time     += wall_t2 - wall_t1;
  OP_kernels[4].transfer += (float)set->size * arg1.size;
  OP_kernels[4].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg0.size;
  OP_kernels[4].transfer += (float)set->size * arg1.map->dim * 4.0f;
}
