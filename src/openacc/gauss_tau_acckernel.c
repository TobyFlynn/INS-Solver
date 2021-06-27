//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void gauss_tau_openacc( const int *edgeNum, const double **fscale, double **tau) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  if(fscale[0][edgeL * 5] > fscale[1][edgeR * 5]) {
    tau[0][edgeL] += 20 * 25 * fscale[0][edgeL * 5];
    tau[1][edgeR] += 20 * 25 * fscale[0][edgeL * 5];
  } else {
    tau[0][edgeL] += 20 * 25 * fscale[1][edgeR * 5];
    tau[1][edgeR] += 20 * 25 * fscale[1][edgeR * 5];
  }
}

// host stub function
void op_par_loop_gauss_tau(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  arg1.idx = 0;
  args[1] = arg1;
  for ( int v=1; v<2; v++ ){
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 15, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 3, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(9);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[9].name      = name;
  OP_kernels[9].count    += 1;

  int  ninds   = 2;
  int  inds[5] = {-1,0,0,1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_tau\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_9
    int part_size = OP_PART_SIZE_9;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map1 = arg1.map_data_d;

    int* data0 = (int*)arg0.data_d;
    double *data1 = (double *)arg1.data_d;
    double *data3 = (double *)arg3.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map1,data0,data1,data3)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map1idx;
        int map2idx;
        map1idx = map1[n + set_size1 * 0];
        map2idx = map1[n + set_size1 * 1];

        const double* arg1_vec[] = {
           &data1[15 * map1idx],
           &data1[15 * map2idx]};
        double* arg3_vec[] = {
           &data3[3 * map1idx],
           &data3[3 * map2idx]};

        gauss_tau_openacc(
          &data0[2 * n],
          arg1_vec,
          arg3_vec);
      }

    }
    OP_kernels[9].transfer  += Plan->transfer;
    OP_kernels[9].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[9].time     += wall_t2 - wall_t1;
}
