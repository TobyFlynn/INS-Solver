//
// auto-generated by op2.py
//

//user function
#include "../kernels/glb_ind_kernelBC.h"

// host stub function
void op_par_loop_glb_ind_kernelBC(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(12);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: glb_ind_kernelBC\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map0idx;
      map0idx = arg0.map_data[n * arg0.map->dim + 0];


      glb_ind_kernelBC(
        &((int*)arg0.data)[1 * map0idx],
        &((int*)arg1.data)[1 * n]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[12].name      = name;
  OP_kernels[12].count    += 1;
  OP_kernels[12].time     += wall_t2 - wall_t1;
  OP_kernels[12].transfer += (float)set->size * arg0.size;
  OP_kernels[12].transfer += (float)set->size * arg1.size;
  OP_kernels[12].transfer += (float)set->size * arg0.map->dim * 4.0f;
}
