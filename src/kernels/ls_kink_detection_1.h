inline void ls_kink_detection_1(const DG_FP **kink_at_nodes, DG_FP *kink) {
  if(kink_at_nodes[0][0] != 0.0 || kink_at_nodes[1][0] != 0.0 || kink_at_nodes[2][0] != 0.0) {
    for(int i = 0; i < DG_NP; i++) {
      kink[i] = 1.0;
    }
  } else {
    for(int i = 0; i < DG_NP; i++) {
      kink[i] = 0.0;
    }
  }
}