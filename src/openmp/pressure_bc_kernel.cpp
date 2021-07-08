//
// auto-generated by op2.py
//

//user function
#include "../kernels/pressure_bc.h"

// host stub function
void op_par_loop_pressure_bc(char const *name, op_set set,
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
  op_arg arg14){

  int nargs = 15;
  op_arg args[15];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(34);
  OP_kernels[34].name      = name;
  OP_kernels[34].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 11;
  int  inds[15] = {-1,-1,-1,-1,0,1,2,3,4,5,6,7,8,9,10};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_bc\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_34
    int part_size = OP_PART_SIZE_34;
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
          int map4idx;
          map4idx = arg4.map_data[n * arg4.map->dim + 0];


          pressure_bc(
            &((int*)arg0.data)[1 * n],
            &((int*)arg1.data)[1 * n],
            (double*)arg2.data,
            (int*)arg3.data,
            &((double*)arg4.data)[15 * map4idx],
            &((double*)arg5.data)[15 * map4idx],
            &((double*)arg6.data)[15 * map4idx],
            &((double*)arg7.data)[15 * map4idx],
            &((double*)arg8.data)[15 * map4idx],
            &((double*)arg9.data)[15 * map4idx],
            &((double*)arg10.data)[15 * map4idx],
            &((double*)arg11.data)[15 * map4idx],
            &((double*)arg12.data)[15 * map4idx],
            &((double*)arg13.data)[15 * map4idx],
            &((double*)arg14.data)[15 * map4idx]);
        }
      }

      block_offset += nblocks;
    }
    OP_kernels[34].transfer  += Plan->transfer;
    OP_kernels[34].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[34].time     += wall_t2 - wall_t1;
}
