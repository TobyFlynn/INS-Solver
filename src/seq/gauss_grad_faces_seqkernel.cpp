//
// auto-generated by op2.py
//

//user function
#include "../kernels/gauss_grad_faces.h"

// host stub function
void op_par_loop_gauss_grad_faces(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5,
  op_arg arg7,
  op_arg arg9,
  op_arg arg11,
  op_arg arg13,
  op_arg arg15,
  op_arg arg17,
  op_arg arg19,
  op_arg arg21,
  op_arg arg23){

  int nargs = 25;
  op_arg args[25];

  args[0] = arg0;
  arg1.idx = 0;
  args[1] = arg1;
  for ( int v=1; v<2; v++ ){
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 105, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 105, "double", OP_READ);
  }

  arg5.idx = 0;
  args[5] = arg5;
  for ( int v=1; v<2; v++ ){
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 105, "double", OP_READ);
  }

  arg7.idx = 0;
  args[7] = arg7;
  for ( int v=1; v<2; v++ ){
    args[7 + v] = op_arg_dat(arg7.dat, v, arg7.map, 105, "double", OP_READ);
  }

  arg9.idx = 0;
  args[9] = arg9;
  for ( int v=1; v<2; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, 105, "double", OP_READ);
  }

  arg11.idx = 0;
  args[11] = arg11;
  for ( int v=1; v<2; v++ ){
    args[11 + v] = op_arg_dat(arg11.dat, v, arg11.map, 105, "double", OP_READ);
  }

  arg13.idx = 0;
  args[13] = arg13;
  for ( int v=1; v<2; v++ ){
    args[13 + v] = op_arg_dat(arg13.dat, v, arg13.map, 105, "double", OP_INC);
  }

  arg15.idx = 0;
  args[15] = arg15;
  for ( int v=1; v<2; v++ ){
    args[15 + v] = op_arg_dat(arg15.dat, v, arg15.map, 105, "double", OP_INC);
  }

  arg17.idx = 0;
  args[17] = arg17;
  for ( int v=1; v<2; v++ ){
    args[17 + v] = op_arg_dat(arg17.dat, v, arg17.map, 105, "double", OP_INC);
  }

  arg19.idx = 0;
  args[19] = arg19;
  for ( int v=1; v<2; v++ ){
    args[19 + v] = op_arg_dat(arg19.dat, v, arg19.map, 105, "double", OP_INC);
  }

  arg21.idx = 0;
  args[21] = arg21;
  for ( int v=1; v<2; v++ ){
    args[21 + v] = op_arg_dat(arg21.dat, v, arg21.map, 105, "double", OP_INC);
  }

  arg23.idx = 0;
  args[23] = arg23;
  for ( int v=1; v<2; v++ ){
    args[23 + v] = op_arg_dat(arg23.dat, v, arg23.map, 105, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(9);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_grad_faces\n");
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
         &((double*)arg1.data)[105 * map1idx],
         &((double*)arg1.data)[105 * map2idx]};
      const double* arg3_vec[] = {
         &((double*)arg3.data)[105 * map1idx],
         &((double*)arg3.data)[105 * map2idx]};
      const double* arg5_vec[] = {
         &((double*)arg5.data)[105 * map1idx],
         &((double*)arg5.data)[105 * map2idx]};
      const double* arg7_vec[] = {
         &((double*)arg7.data)[105 * map1idx],
         &((double*)arg7.data)[105 * map2idx]};
      const double* arg9_vec[] = {
         &((double*)arg9.data)[105 * map1idx],
         &((double*)arg9.data)[105 * map2idx]};
      const double* arg11_vec[] = {
         &((double*)arg11.data)[105 * map1idx],
         &((double*)arg11.data)[105 * map2idx]};
      double* arg13_vec[] = {
         &((double*)arg13.data)[105 * map1idx],
         &((double*)arg13.data)[105 * map2idx]};
      double* arg15_vec[] = {
         &((double*)arg15.data)[105 * map1idx],
         &((double*)arg15.data)[105 * map2idx]};
      double* arg17_vec[] = {
         &((double*)arg17.data)[105 * map1idx],
         &((double*)arg17.data)[105 * map2idx]};
      double* arg19_vec[] = {
         &((double*)arg19.data)[105 * map1idx],
         &((double*)arg19.data)[105 * map2idx]};
      double* arg21_vec[] = {
         &((double*)arg21.data)[105 * map1idx],
         &((double*)arg21.data)[105 * map2idx]};
      double* arg23_vec[] = {
         &((double*)arg23.data)[105 * map1idx],
         &((double*)arg23.data)[105 * map2idx]};

      gauss_grad_faces(
        &((int*)arg0.data)[2 * n],
        arg1_vec,
        arg3_vec,
        arg5_vec,
        arg7_vec,
        arg9_vec,
        arg11_vec,
        arg13_vec,
        arg15_vec,
        arg17_vec,
        arg19_vec,
        arg21_vec,
        arg23_vec);
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
  OP_kernels[9].transfer += (float)set->size * arg1.size;
  OP_kernels[9].transfer += (float)set->size * arg3.size;
  OP_kernels[9].transfer += (float)set->size * arg5.size;
  OP_kernels[9].transfer += (float)set->size * arg7.size;
  OP_kernels[9].transfer += (float)set->size * arg9.size;
  OP_kernels[9].transfer += (float)set->size * arg11.size;
  OP_kernels[9].transfer += (float)set->size * arg13.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg15.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg17.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg19.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg21.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg23.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg0.size;
  OP_kernels[9].transfer += (float)set->size * arg1.map->dim * 4.0f;
}
