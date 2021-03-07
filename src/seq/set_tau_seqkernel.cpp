//
// auto-generated by op2.py
//

//user function
#include "../kernels/set_tau.h"

// host stub function
void op_par_loop_set_tau(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5,
  op_arg arg7,
  op_arg arg9){

  int nargs = 11;
  op_arg args[11];

  args[0] = arg0;
  arg1.idx = 0;
  args[1] = arg1;
  for ( int v=1; v<2; v++ ){
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 3, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 3, "double", OP_READ);
  }

  arg5.idx = 0;
  args[5] = arg5;
  for ( int v=1; v<2; v++ ){
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 15, "double", OP_READ);
  }

  arg7.idx = 0;
  args[7] = arg7;
  for ( int v=1; v<2; v++ ){
    args[7 + v] = op_arg_dat(arg7.dat, v, arg7.map, 15, "double", OP_READ);
  }

  arg9.idx = 0;
  args[9] = arg9;
  for ( int v=1; v<2; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, 15, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(19);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: set_tau\n");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map1idx;
      int map2idx;
      map1idx = arg1.map_data[n * arg1.map->dim + 0];
      map2idx = arg1.map_data[n * arg1.map->dim + 1];

      const double* arg1_vec[] = {
         &((double*)arg1.data)[3 * map1idx],
         &((double*)arg1.data)[3 * map2idx]};
      const double* arg3_vec[] = {
         &((double*)arg3.data)[3 * map1idx],
         &((double*)arg3.data)[3 * map2idx]};
      const double* arg5_vec[] = {
         &((double*)arg5.data)[15 * map1idx],
         &((double*)arg5.data)[15 * map2idx]};
      const double* arg7_vec[] = {
         &((double*)arg7.data)[15 * map1idx],
         &((double*)arg7.data)[15 * map2idx]};
      double* arg9_vec[] = {
         &((double*)arg9.data)[15 * map1idx],
         &((double*)arg9.data)[15 * map2idx]};

      set_tau(
        &((int*)arg0.data)[2 * n],
        arg1_vec,
        arg3_vec,
        arg5_vec,
        arg7_vec,
        arg9_vec);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[19].name      = name;
  OP_kernels[19].count    += 1;
  OP_kernels[19].time     += wall_t2 - wall_t1;
  OP_kernels[19].transfer += (float)set->size * arg1.size;
  OP_kernels[19].transfer += (float)set->size * arg3.size;
  OP_kernels[19].transfer += (float)set->size * arg5.size;
  OP_kernels[19].transfer += (float)set->size * arg7.size;
  OP_kernels[19].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[19].transfer += (float)set->size * arg0.size;
  OP_kernels[19].transfer += (float)set->size * arg1.map->dim * 4.0f;
}
