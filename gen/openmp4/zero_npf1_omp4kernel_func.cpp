//
// auto-generated by op2.py
//

void zero_npf1_omp4_kernel(
  double *data0,
  int dat0size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    double *npf0 = &data0[12*n_op];

    //inline function
    
    for(int i = 0; i < 3 * 4; i++) {
      npf0[i] = 0.0;
    }
    //end inline func
  }

}
