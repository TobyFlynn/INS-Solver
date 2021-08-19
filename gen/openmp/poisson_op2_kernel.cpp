//
// auto-generated by op2.py
//

//user function
#include "../kernels/poisson_op2.h"

// host stub function
void op_par_loop_poisson_op2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg17,
  op_arg arg18,
  op_arg arg19){

  int nargs = 20;
  op_arg args[20];

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
  args[10] = arg10;
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;
  args[16] = arg16;
  args[17] = arg17;
  args[18] = arg18;
  args[19] = arg19;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(17);
  OP_kernels[17].name      = name;
  OP_kernels[17].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 5;
  int  inds[20] = {-1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,2,2,3,3,4,4,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op2\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_17
    int part_size = OP_PART_SIZE_17;
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
          int map8idx;
          int map9idx;
          map8idx = arg8.map_data[n * arg8.map->dim + 0];
          map9idx = arg8.map_data[n * arg8.map->dim + 1];


          poisson_op2(
            &((int*)arg0.data)[2 * n],
            &((bool*)arg1.data)[1 * n],
            &((double*)arg2.data)[60 * n],
            &((double*)arg3.data)[60 * n],
            &((double*)arg4.data)[60 * n],
            &((double*)arg5.data)[60 * n],
            &((double*)arg6.data)[60 * n],
            &((double*)arg7.data)[60 * n],
            &((double*)arg8.data)[18 * map8idx],
            &((double*)arg8.data)[18 * map9idx],
            &((double*)arg10.data)[1 * map8idx],
            &((double*)arg10.data)[1 * map9idx],
            &((double*)arg12.data)[18 * map8idx],
            &((double*)arg12.data)[18 * map9idx],
            &((double*)arg14.data)[10 * map8idx],
            &((double*)arg14.data)[10 * map9idx],
            &((double*)arg16.data)[100 * map8idx],
            &((double*)arg16.data)[100 * map9idx],
            &((double*)arg18.data)[100 * n],
            &((double*)arg19.data)[100 * n]);
        }
      }

      block_offset += nblocks;
    }
    OP_kernels[17].transfer  += Plan->transfer;
    OP_kernels[17].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[17].time     += wall_t2 - wall_t1;
}
