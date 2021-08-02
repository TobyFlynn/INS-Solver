//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void sigma_flux_openacc( const int *edgeNum, const bool *rev, const double **sJ,
                       const double **nx, const double **ny, const double **s,
                       double **sigFx, double **sigFy) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exIndL = edgeL * 6;
  int exIndR = edgeR * 6;

  for(int i = 0; i < 6; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + 6 - i - 1;
    } else {
      rInd = exIndR + i;
    }
    double flux = (s[0][lInd] + s[1][rInd]) / 2.0;
    sigFx[0][lInd] += gaussW_g[i] * sJ[0][lInd] * nx[0][lInd] * flux;
    sigFy[0][lInd] += gaussW_g[i] * sJ[0][lInd] * ny[0][lInd] * flux;
  }

  for(int i = 0; i < 6; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 6 - i - 1;
    } else {
      lInd = exIndL + i;
    }
    double flux = (s[0][lInd] + s[1][rInd]) / 2.0;
    sigFx[1][rInd] += gaussW_g[i] * sJ[1][rInd] * nx[1][rInd] * flux;
    sigFy[1][rInd] += gaussW_g[i] * sJ[1][rInd] * ny[1][rInd] * flux;
  }
}

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
  op_timing_realloc(58);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[58].name      = name;
  OP_kernels[58].count    += 1;

  int  ninds   = 6;
  int  inds[14] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: sigma_flux\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_58
    int part_size = OP_PART_SIZE_58;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map2 = arg2.map_data_d;

    int* data0 = (int*)arg0.data_d;
    bool* data1 = (bool*)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data4 = (double *)arg4.data_d;
    double *data6 = (double *)arg6.data_d;
    double *data8 = (double *)arg8.data_d;
    double *data10 = (double *)arg10.data_d;
    double *data12 = (double *)arg12.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data4,data6,data8,data10,data12)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        int map3idx;
        map2idx = map2[n + set_size1 * 0];
        map3idx = map2[n + set_size1 * 1];

        const double* arg2_vec[] = {
           &data2[18 * map2idx],
           &data2[18 * map3idx]};
        const double* arg4_vec[] = {
           &data4[18 * map2idx],
           &data4[18 * map3idx]};
        const double* arg6_vec[] = {
           &data6[18 * map2idx],
           &data6[18 * map3idx]};
        const double* arg8_vec[] = {
           &data8[18 * map2idx],
           &data8[18 * map3idx]};
        double* arg10_vec[] = {
           &data10[18 * map2idx],
           &data10[18 * map3idx]};
        double* arg12_vec[] = {
           &data12[18 * map2idx],
           &data12[18 * map3idx]};

        sigma_flux_openacc(
          &data0[2 * n],
          &data1[1 * n],
          arg2_vec,
          arg4_vec,
          arg6_vec,
          arg8_vec,
          arg10_vec,
          arg12_vec);
      }

    }
    OP_kernels[58].transfer  += Plan->transfer;
    OP_kernels[58].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[58].time     += wall_t2 - wall_t1;
}
