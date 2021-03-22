inline void min_max(double *min, double *max, const double *arr) {
  for(int i = 0; i < 15; i++) {
    if(arr[i] > *max) {
      *max = arr[i];
    }
    if(arr[i] < *min) {
      *min = arr[i];
    }
  }
}
