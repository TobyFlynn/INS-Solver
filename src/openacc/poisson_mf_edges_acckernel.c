//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_mf_edges_openacc( const int *edgeNum, const bool *rev,
                             const double **sJ, const double **nx,
                             const double **ny, const double **tau,
                             const double **u, const double **dudx,
                             const double **dudy, double **fluxX,
                             double **fluxY, double **flux) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exIndL = 0;
  if(edgeL == 1) exIndL = 7;
  else if(edgeL == 2) exIndL = 2 * 7;

  int exIndR = 0;
  if(edgeR == 1) exIndR = 7;
  else if(edgeR == 2) exIndR = 2 * 7;

  for(int i = 0; i < 7; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + 7 - i - 1;
    } else {
      rInd = exIndR + i;
    }

    double tmp = (u[0][lInd] + u[1][rInd]) / 2.0;
    tmp *= gaussW_g[i] * sJ[0][lInd];
    fluxX[0][lInd] += nx[0][lInd] * tmp;
    fluxY[0][lInd] += ny[0][lInd] * tmp;
    tmp = nx[0][lInd] * ((dudx[0][lInd] + dudx[1][rInd]) / 2.0);
    tmp += ny[0][lInd] * ((dudy[0][lInd] + dudy[1][rInd]) / 2.0);
    tmp -= tau[0][edgeL] * (u[0][lInd] - u[1][rInd]) / 2.0;
    tmp *= gaussW_g[i] * sJ[0][lInd];
    flux[0][lInd] += tmp;
  }

  for(int i = 0; i < 7; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 7 - i - 1;
    } else {
      lInd = exIndL + i;
    }

    double tmp = (u[0][lInd] + u[1][rInd]) / 2.0;
    tmp *= gaussW_g[i] * sJ[1][rInd];
    fluxX[1][rInd] += nx[1][rInd] * tmp;
    fluxY[1][rInd] += ny[1][rInd] * tmp;
    tmp = nx[1][rInd] * ((dudx[0][lInd] + dudx[1][rInd]) / 2.0);
    tmp += ny[1][rInd] * ((dudy[0][lInd] + dudy[1][rInd]) / 2.0);
    tmp -= tau[1][edgeR] * (u[1][rInd] - u[0][lInd]) / 2.0;
    tmp *= gaussW_g[i] * sJ[1][rInd];
    flux[1][rInd] += tmp;
  }
}

// host stub function
void op_par_loop_poisson_mf_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6,
  op_arg arg8,
  op_arg arg10,
  op_arg arg12,
  op_arg arg14,
  op_arg arg16,
  op_arg arg18,
  op_arg arg20){

  int nargs = 22;
  op_arg args[22];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 21, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 21, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 21, "double", OP_READ);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 3, "double", OP_READ);
  }

  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 21, "double", OP_READ);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 21, "double", OP_READ);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<2; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 21, "double", OP_READ);
  }

  arg16.idx = 0;
  args[16] = arg16;
  for ( int v=1; v<2; v++ ){
    args[16 + v] = op_arg_dat(arg16.dat, v, arg16.map, 21, "double", OP_INC);
  }

  arg18.idx = 0;
  args[18] = arg18;
  for ( int v=1; v<2; v++ ){
    args[18 + v] = op_arg_dat(arg18.dat, v, arg18.map, 21, "double", OP_INC);
  }

  arg20.idx = 0;
  args[20] = arg20;
  for ( int v=1; v<2; v++ ){
    args[20 + v] = op_arg_dat(arg20.dat, v, arg20.map, 21, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(20);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[20].name      = name;
  OP_kernels[20].count    += 1;

  int  ninds   = 10;
  int  inds[22] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_mf_edges\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_20
    int part_size = OP_PART_SIZE_20;
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
    double *data14 = (double *)arg14.data_d;
    double *data16 = (double *)arg16.data_d;
    double *data18 = (double *)arg18.data_d;
    double *data20 = (double *)arg20.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data4,data6,data8,data10,data12,data14,data16,data18,data20)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        int map3idx;
        map2idx = map2[n + set_size1 * 0];
        map3idx = map2[n + set_size1 * 1];

        const double* arg2_vec[] = {
           &data2[21 * map2idx],
           &data2[21 * map3idx]};
        const double* arg4_vec[] = {
           &data4[21 * map2idx],
           &data4[21 * map3idx]};
        const double* arg6_vec[] = {
           &data6[21 * map2idx],
           &data6[21 * map3idx]};
        const double* arg8_vec[] = {
           &data8[3 * map2idx],
           &data8[3 * map3idx]};
        const double* arg10_vec[] = {
           &data10[21 * map2idx],
           &data10[21 * map3idx]};
        const double* arg12_vec[] = {
           &data12[21 * map2idx],
           &data12[21 * map3idx]};
        const double* arg14_vec[] = {
           &data14[21 * map2idx],
           &data14[21 * map3idx]};
        double* arg16_vec[] = {
           &data16[21 * map2idx],
           &data16[21 * map3idx]};
        double* arg18_vec[] = {
           &data18[21 * map2idx],
           &data18[21 * map3idx]};
        double* arg20_vec[] = {
           &data20[21 * map2idx],
           &data20[21 * map3idx]};

        poisson_mf_edges_openacc(
          &data0[2 * n],
          &data1[1 * n],
          arg2_vec,
          arg4_vec,
          arg6_vec,
          arg8_vec,
          arg10_vec,
          arg12_vec,
          arg14_vec,
          arg16_vec,
          arg18_vec,
          arg20_vec);
      }

    }
    OP_kernels[20].transfer  += Plan->transfer;
    OP_kernels[20].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[20].time     += wall_t2 - wall_t1;
}
