//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void glb_ind_kernel_openacc( const int **glb, int *glbL, int *glbR) {
  glbL[0] = glb[0][0];
  glbR[0] = glb[1][0];
}

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
  op_timing_realloc(23);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[23].name      = name;
  OP_kernels[23].count    += 1;

  int  ninds   = 1;
  int  inds[4] = {0,0,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: glb_ind_kernel\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_23
    int part_size = OP_PART_SIZE_23;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map0 = arg0.map_data_d;

    int* data2 = (int*)arg2.data_d;
    int* data3 = (int*)arg3.data_d;
    int *data0 = (int *)arg0.data_d;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;
    int set_size1 = set->size + set->exec_size;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      #pragma acc parallel loop independent deviceptr(col_reord,map0,data2,data3,data0)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map0idx;
        int map1idx;
        map0idx = map0[n + set_size1 * 0];
        map1idx = map0[n + set_size1 * 1];

        const int* arg0_vec[] = {
           &data0[1 * map0idx],
           &data0[1 * map1idx]};

        glb_ind_kernel_openacc(
          arg0_vec,
          &data2[1 * n],
          &data3[1 * n]);
      }

    }
    OP_kernels[23].transfer  += Plan->transfer;
    OP_kernels[23].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[23].time     += wall_t2 - wall_t1;
}
