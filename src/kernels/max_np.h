inline void max_np(const DG_FP *data, DG_FP *max) {
  DG_FP local_max = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    if(local_max < data[i])
      local_max = data[i];
  }
  if(*max < local_max)
    *max = local_max;
}
