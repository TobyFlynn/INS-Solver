//
// auto-generated by op2.py
//

//user function
#include "../kernels/ls_advec_edges.h"

// host stub function
void op_par_loop_ls_advec_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 15, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 15, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(46);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_advec_edges\n");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map2idx;
      int map3idx;
      map2idx = arg2.map_data[n * arg2.map->dim + 0];
      map3idx = arg2.map_data[n * arg2.map->dim + 1];

      const double* arg2_vec[] = {
         &((double*)arg2.data)[15 * map2idx],
         &((double*)arg2.data)[15 * map3idx]};
      double* arg4_vec[] = {
         &((double*)arg4.data)[15 * map2idx],
         &((double*)arg4.data)[15 * map3idx]};

      ls_advec_edges(
        &((int*)arg0.data)[2 * n],
        &((bool*)arg1.data)[1 * n],
        arg2_vec,
        arg4_vec);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[46].name      = name;
  OP_kernels[46].count    += 1;
  OP_kernels[46].time     += wall_t2 - wall_t1;
  OP_kernels[46].transfer += (float)set->size * arg2.size;
  OP_kernels[46].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[46].transfer += (float)set->size * arg0.size;
  OP_kernels[46].transfer += (float)set->size * arg1.size;
  OP_kernels[46].transfer += (float)set->size * arg2.map->dim * 4.0f;
}
