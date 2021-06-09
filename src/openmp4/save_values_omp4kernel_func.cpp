//
// auto-generated by op2.py
//

void save_values_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *v_vals = &data0[15*n_op];
    double *c_vals = &data1[16*n_op];

    //inline function
    
    c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[5]) / 3.0;
    c_vals[1]  = (v_vals[1] + v_vals[5] + v_vals[6]) / 3.0;
    c_vals[2]  = (v_vals[1] + v_vals[2] + v_vals[6]) / 3.0;
    c_vals[3]  = (v_vals[2] + v_vals[6] + v_vals[7]) / 3.0;
    c_vals[4]  = (v_vals[2] + v_vals[3] + v_vals[7]) / 3.0;
    c_vals[5]  = (v_vals[3] + v_vals[7] + v_vals[8]) / 3.0;
    c_vals[6]  = (v_vals[3] + v_vals[4] + v_vals[8]) / 3.0;
    c_vals[7]  = (v_vals[5] + v_vals[6] + v_vals[9]) / 3.0;
    c_vals[8]  = (v_vals[6] + v_vals[9] + v_vals[10]) / 3.0;
    c_vals[9]  = (v_vals[6] + v_vals[7] + v_vals[10]) / 3.0;
    c_vals[10] = (v_vals[7] + v_vals[10] + v_vals[11]) / 3.0;
    c_vals[11] = (v_vals[7] + v_vals[8] + v_vals[11]) / 3.0;
    c_vals[12] = (v_vals[9] + v_vals[10] + v_vals[12]) / 3.0;
    c_vals[13] = (v_vals[10] + v_vals[12] + v_vals[13]) / 3.0;
    c_vals[14] = (v_vals[10] + v_vals[11] + v_vals[13]) / 3.0;
    c_vals[15] = (v_vals[12] + v_vals[13] + v_vals[14]) / 3.0;
    //end inline func
  }

}