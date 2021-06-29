//
// auto-generated by op2.py
//

//user function
#include "../kernels/pressure_solve_1.h"

// host stub function
void op_par_loop_pressure_solve_1(char const *name, op_set set,
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
  op_arg arg19,
  op_arg arg20,
  op_arg arg21,
  op_arg arg22,
  op_arg arg23,
  op_arg arg24,
  op_arg arg25,
  op_arg arg26,
  op_arg arg27,
  op_arg arg28,
  op_arg arg29,
  op_arg arg30,
  op_arg arg31,
  op_arg arg32,
  op_arg arg33){

  int nargs = 34;
  op_arg args[34];

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
  args[20] = arg20;
  args[21] = arg21;
  args[22] = arg22;
  args[23] = arg23;
  args[24] = arg24;
  args[25] = arg25;
  args[26] = arg26;
  args[27] = arg27;
  args[28] = arg28;
  args[29] = arg29;
  args[30] = arg30;
  args[31] = arg31;
  args[32] = arg32;
  args[33] = arg33;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(28);
  OP_kernels[28].name      = name;
  OP_kernels[28].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 16;
  int  inds[34] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_solve_1\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_28
    int part_size = OP_PART_SIZE_28;
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
          int map2idx;
          int map3idx;
          map2idx = arg2.map_data[n * arg2.map->dim + 0];
          map3idx = arg2.map_data[n * arg2.map->dim + 1];


          pressure_solve_1(
            &((int*)arg0.data)[2 * n],
            &((bool*)arg1.data)[1 * n],
            &((double*)arg2.data)[105 * map2idx],
            &((double*)arg2.data)[105 * map3idx],
            &((double*)arg4.data)[105 * map2idx],
            &((double*)arg4.data)[105 * map3idx],
            &((double*)arg6.data)[105 * map2idx],
            &((double*)arg6.data)[105 * map3idx],
            &((double*)arg8.data)[105 * map2idx],
            &((double*)arg8.data)[105 * map3idx],
            &((double*)arg10.data)[105 * map2idx],
            &((double*)arg10.data)[105 * map3idx],
            &((double*)arg12.data)[105 * map2idx],
            &((double*)arg12.data)[105 * map3idx],
            &((double*)arg14.data)[105 * map2idx],
            &((double*)arg14.data)[105 * map3idx],
            &((double*)arg16.data)[105 * map2idx],
            &((double*)arg16.data)[105 * map3idx],
            &((double*)arg18.data)[105 * map2idx],
            &((double*)arg18.data)[105 * map3idx],
            &((double*)arg20.data)[21 * map2idx],
            &((double*)arg20.data)[21 * map3idx],
            &((double*)arg22.data)[1 * map2idx],
            &((double*)arg22.data)[1 * map3idx],
            &((double*)arg24.data)[3 * map2idx],
            &((double*)arg24.data)[3 * map3idx],
            &((double*)arg26.data)[21 * map2idx],
            &((double*)arg26.data)[21 * map3idx],
            &((double*)arg28.data)[15 * map2idx],
            &((double*)arg28.data)[15 * map3idx],
            &((double*)arg30.data)[15 * map2idx],
            &((double*)arg30.data)[15 * map3idx],
            &((double*)arg32.data)[15 * map2idx],
            &((double*)arg32.data)[15 * map3idx]);
        }
      }

      block_offset += nblocks;
    }
    OP_kernels[28].transfer  += Plan->transfer;
    OP_kernels[28].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[28].time     += wall_t2 - wall_t1;
}
