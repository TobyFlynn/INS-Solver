//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void pressure_solve_1_openacc( const int *edgeNum, const bool *rev,
                             const double **mD0, const double **mD1,
                             const double **mD2, const double **pD0,
                             const double **pD1, const double **pD2,
                             const double **gVP0, const double **gVP1,
                             const double **gVP2, const double **sJ,
                             const double **h, const double **tau, const double **rho,
                             const double **u, double *rhsL, double *rhsR) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;


  const double *mDL, *mDR, *pDL, *pDR, *gVML, *gVMR, *gVPL, *gVPR;
  if(edgeL == 0) {
    mDL  = mD0[0];
    pDL  = pD0[0];
    gVML = gFInterp0_g;
    gVPL = gVP0[0];
  } else if(edgeL == 1) {
    mDL  = mD1[0];
    pDL  = pD1[0];
    gVML = gFInterp1_g;
    gVPL = gVP1[0];
  } else {
    mDL  = mD2[0];
    pDL  = pD2[0];
    gVML = gFInterp2_g;
    gVPL = gVP2[0];
  }

  if(edgeR == 0) {
    mDR  = mD0[1];
    pDR  = pD0[1];
    gVMR = gFInterp0_g;
    gVPR = gVP0[1];
  } else if(edgeR == 1) {
    mDR  = mD1[1];
    pDR  = pD1[1];
    gVMR = gFInterp1_g;
    gVPR = gVP1[1];
  } else {
    mDR  = mD2[1];
    pDR  = pD2[1];
    gVMR = gFInterp2_g;
    gVPR = gVP2[1];
  }

  double op1L[15 * 15];
  double op1R[15 * 15];
  double op2L[15 * 15];
  double op2R[15 * 15];



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      op1L[c_ind] = 0.0;
      op1R[c_ind] = 0.0;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);


        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;

        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * (1.0 / rho[0][factors_indL]) * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * (1.0 / rho[1][factors_indR]) * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * (1.0 / rho[0][factors_indL]) * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * (1.0 / rho[1][factors_indR]) * pDR[b_ind];
      }
    }
  }



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);

        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;










        op1L[c_ind] += -(1.0 / rho[0][factors_indL]) * mDL[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVML[b_ind];
        op1R[c_ind] += -(1.0 / rho[1][factors_indR]) * mDR[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVMR[b_ind];

        op2L[c_ind] += (1.0 / rho[0][factors_indL]) * mDL[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVPL[b_ind];
        op2R[c_ind] += (1.0 / rho[1][factors_indR]) * mDR[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVPR[b_ind];
      }
    }
  }

  double tauL[7];
  double tauR[7];
  for(int i = 0; i < 7; i++) {
    int indL = edgeL * 7 + i;
    int indR;
    if(reverse)
      indR = edgeR * 7 + 6 - i;
    else
      indR = edgeR * 7 + i;
    tauL[i] = 10 * 0.5 * 5 * 6 * fmax(*(h[0]) / rho[0][indL], *(h[1]) / rho[1][indR]);
  }
  for(int i = 0; i < 7; i++) {
    int indL;
    int indR = edgeR * 7 + i;
    if(reverse)
      indL = edgeL * 7 + 6 - i;
    else
      indL = edgeL * 7 + i;
    tauR[i] = 10 * 0.5 * 5 * 6 * fmax(*(h[0]) / rho[0][indL], *(h[1]) / rho[1][indR]);
  }



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);

        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;










        op1L[c_ind] += gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * tauR[k] * gVPR[b_ind];
      }
    }
  }


  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int op_ind = i * 15 + j;
      rhsL[i] += op1L[op_ind] * u[0][j] + op2L[op_ind] * u[1][j];
      rhsR[i] += op1R[op_ind] * u[1][j] + op2R[op_ind] * u[0][j];
    }
  }
}

// host stub function
void op_par_loop_pressure_solve_1(char const *name, op_set set,
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
  op_arg arg20,
  op_arg arg22,
  op_arg arg24,
  op_arg arg26,
  op_arg arg28,
  op_arg arg30,
  op_arg arg31){

  int nargs = 32;
  op_arg args[32];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 105, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 105, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 105, "double", OP_READ);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 105, "double", OP_READ);
  }

  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 105, "double", OP_READ);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 105, "double", OP_READ);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<2; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 105, "double", OP_READ);
  }

  arg16.idx = 0;
  args[16] = arg16;
  for ( int v=1; v<2; v++ ){
    args[16 + v] = op_arg_dat(arg16.dat, v, arg16.map, 105, "double", OP_READ);
  }

  arg18.idx = 0;
  args[18] = arg18;
  for ( int v=1; v<2; v++ ){
    args[18 + v] = op_arg_dat(arg18.dat, v, arg18.map, 105, "double", OP_READ);
  }

  arg20.idx = 0;
  args[20] = arg20;
  for ( int v=1; v<2; v++ ){
    args[20 + v] = op_arg_dat(arg20.dat, v, arg20.map, 21, "double", OP_READ);
  }

  arg22.idx = 0;
  args[22] = arg22;
  for ( int v=1; v<2; v++ ){
    args[22 + v] = op_arg_dat(arg22.dat, v, arg22.map, 1, "double", OP_READ);
  }

  arg24.idx = 0;
  args[24] = arg24;
  for ( int v=1; v<2; v++ ){
    args[24 + v] = op_arg_dat(arg24.dat, v, arg24.map, 3, "double", OP_READ);
  }

  arg26.idx = 0;
  args[26] = arg26;
  for ( int v=1; v<2; v++ ){
    args[26 + v] = op_arg_dat(arg26.dat, v, arg26.map, 21, "double", OP_READ);
  }

  arg28.idx = 0;
  args[28] = arg28;
  for ( int v=1; v<2; v++ ){
    args[28 + v] = op_arg_dat(arg28.dat, v, arg28.map, 15, "double", OP_READ);
  }

  args[30] = arg30;
  args[31] = arg31;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(28);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[28].name      = name;
  OP_kernels[28].count    += 1;

  int  ninds   = 15;
  int  inds[32] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_solve_1\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_28
    int part_size = OP_PART_SIZE_28;
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
    double *data22 = (double *)arg22.data_d;
    double *data24 = (double *)arg24.data_d;
    double *data26 = (double *)arg26.data_d;
    double *data28 = (double *)arg28.data_d;
    double *data30 = (double *)arg30.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data4,data6,data8,data10,data12,data14,data16,data18,data20,data22,data24,data26,data28,data30)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        int map3idx;
        map2idx = map2[n + set_size1 * 0];
        map3idx = map2[n + set_size1 * 1];

        const double* arg2_vec[] = {
           &data2[105 * map2idx],
           &data2[105 * map3idx]};
        const double* arg4_vec[] = {
           &data4[105 * map2idx],
           &data4[105 * map3idx]};
        const double* arg6_vec[] = {
           &data6[105 * map2idx],
           &data6[105 * map3idx]};
        const double* arg8_vec[] = {
           &data8[105 * map2idx],
           &data8[105 * map3idx]};
        const double* arg10_vec[] = {
           &data10[105 * map2idx],
           &data10[105 * map3idx]};
        const double* arg12_vec[] = {
           &data12[105 * map2idx],
           &data12[105 * map3idx]};
        const double* arg14_vec[] = {
           &data14[105 * map2idx],
           &data14[105 * map3idx]};
        const double* arg16_vec[] = {
           &data16[105 * map2idx],
           &data16[105 * map3idx]};
        const double* arg18_vec[] = {
           &data18[105 * map2idx],
           &data18[105 * map3idx]};
        const double* arg20_vec[] = {
           &data20[21 * map2idx],
           &data20[21 * map3idx]};
        const double* arg22_vec[] = {
           &data22[1 * map2idx],
           &data22[1 * map3idx]};
        const double* arg24_vec[] = {
           &data24[3 * map2idx],
           &data24[3 * map3idx]};
        const double* arg26_vec[] = {
           &data26[21 * map2idx],
           &data26[21 * map3idx]};
        const double* arg28_vec[] = {
           &data28[15 * map2idx],
           &data28[15 * map3idx]};

        pressure_solve_1_openacc(
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
          arg20_vec,
          arg22_vec,
          arg24_vec,
          arg26_vec,
          arg28_vec,
          &data30[15 * map2idx],
          &data30[15 * map3idx]);
      }

    }
    OP_kernels[28].transfer  += Plan->transfer;
    OP_kernels[28].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[28].time     += wall_t2 - wall_t1;
}
