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
  op_timing_realloc(22);
  OP_kernels[22].name      = name;
  OP_kernels[22].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 1;
  int  inds[4] = {0,0,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: glb_ind_kernel\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_22
    int part_size = OP_PART_SIZE_22;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    op_plan *Plan = op_plan_get_stage_upload(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL,0);

    // execute plan
    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all(nargs, args);
      }
      int nblocks = Plan->ncolblk[col];

      #pragma omp parallel for
      for ( int blockIdx=0; blockIdx<nblocks; blockIdx++ ){
        int blockId  = Plan->blkmap[blockIdx + block_offset];
        int nelem    = Plan->nelems[blockId];
        int offset_b = Plan->offset[blockId];
        for ( int n=offset_b; n<offset_b+nelem; n++ ){
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

      block_offset += nblocks;
    }
    OP_kernels[22].transfer  += Plan->transfer;
    OP_kernels[22].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[22].time     += wall_t2 - wall_t1;
}
