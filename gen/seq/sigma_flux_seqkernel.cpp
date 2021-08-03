//
// auto-generated by op2.py
//

//user function
#include "../kernels/sigma_flux.h"

// host stub function
void op_par_loop_sigma_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6,
  op_arg arg8,
  op_arg arg10,
  op_arg arg12){

  int nargs = 14;
  op_arg args[14];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 18, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 18, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 18, "double", OP_READ);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 18, "double", OP_READ);
  }

  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 18, "double", OP_INC);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 18, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(60);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: sigma_flux\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map2idx;
      int map3idx;
      map2idx = arg2.map_data[n * arg2.map->dim + 0];
      map3idx = arg2.map_data[n * arg2.map->dim + 1];

      const double* arg2_vec[] = {
         &((double*)arg2.data)[18 * map2idx],
         &((double*)arg2.data)[18 * map3idx]};
      const double* arg4_vec[] = {
         &((double*)arg4.data)[18 * map2idx],
         &((double*)arg4.data)[18 * map3idx]};
      const double* arg6_vec[] = {
         &((double*)arg6.data)[18 * map2idx],
         &((double*)arg6.data)[18 * map3idx]};
      const double* arg8_vec[] = {
         &((double*)arg8.data)[18 * map2idx],
         &((double*)arg8.data)[18 * map3idx]};
      double* arg10_vec[] = {
         &((double*)arg10.data)[18 * map2idx],
         &((double*)arg10.data)[18 * map3idx]};
      double* arg12_vec[] = {
         &((double*)arg12.data)[18 * map2idx],
         &((double*)arg12.data)[18 * map3idx]};

      sigma_flux(
        &((int*)arg0.data)[2 * n],
        &((bool*)arg1.data)[1 * n],
        arg2_vec,
        arg4_vec,
        arg6_vec,
        arg8_vec,
        arg10_vec,
        arg12_vec);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[60].name      = name;
  OP_kernels[60].count    += 1;
  OP_kernels[60].time     += wall_t2 - wall_t1;
  OP_kernels[60].transfer += (float)set->size * arg2.size;
  OP_kernels[60].transfer += (float)set->size * arg4.size;
  OP_kernels[60].transfer += (float)set->size * arg6.size;
  OP_kernels[60].transfer += (float)set->size * arg8.size;
  OP_kernels[60].transfer += (float)set->size * arg10.size * 2.0f;
  OP_kernels[60].transfer += (float)set->size * arg12.size * 2.0f;
  OP_kernels[60].transfer += (float)set->size * arg0.size;
  OP_kernels[60].transfer += (float)set->size * arg1.size;
  OP_kernels[60].transfer += (float)set->size * arg2.map->dim * 4.0f;
}
