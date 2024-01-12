inline void modal_shock_detector_3d_0(const DG_FP *max_vis, const DG_FP *modes,
                                      DG_FP **vis_nodes, DG_FP **count_nodes) {
  // Compute 1D equivalent of modes
  DG_FP q[DG_ORDER + 1];
  int counts[DG_ORDER + 1];
  for(int i = 0; i < DG_ORDER + 1; i++) {
    q[i] = 0.0;
    counts[i] = 0;
  }

  int ind = 0;
  for(int i = 0; i < DG_ORDER + 1; i++) {
    for(int j = 0; j < DG_ORDER - i + 1; j++) {
      for(int k = 0; k < DG_ORDER - i - j + 1; k++) {
        int order = i + j + k;
        q[order] += modes[ind] * modes[ind];
        counts[order]++;
        ind++;
      }
    }
  }

  for(int i = 0; i < DG_ORDER + 1; i++) {
    // q[i] = sqrt(q[i] / (DG_FP)counts[i]);
    q[i] = sqrt(q[i]);
  }

  // Skyline pessimization
  DG_FP q_bar[DG_ORDER + 1];
  q_bar[0] = q[0];
  for(int i = 1; i < DG_ORDER + 1; i++) {
    q_bar[i] = q[i];
    int order = i < DG_ORDER - 1 ? i : DG_ORDER - 1;
    for(int j = order; j < DG_ORDER + 1; j++) {
      if(q_bar[i] < q[j]) q_bar[i] = q[j];
    }
  }

  // Least squares fit
  DG_FP sum_0 = 0.0;
  DG_FP sum_1 = 0.0;
  DG_FP sum_2 = 0.0;
  DG_FP sum_3 = 0.0;
  for(int i = 1; i < DG_ORDER + 1; i++) {
    const DG_FP log_q = log10(q_bar[i]);
    const DG_FP log_x = log10((DG_FP)i + 1.0);
    sum_0 += log_q * log_x;
    sum_1 += log_q;
    sum_2 += log_x;
    sum_3 += log_x * log_x;
  }
  const DG_FP n = (DG_FP)(DG_ORDER);
  const DG_FP s = -(n * sum_0 - sum_1 * sum_2) / (n * sum_3 - sum_2 * sum_2);

  DG_FP art_vis = 0.0;
  const DG_FP PI = 3.141592653589793238463;
  if(s < 1.0) art_vis = *max_vis;
  else if(s <= 3.0) art_vis = *max_vis * 0.5 * (1 + sin(-PI * (s - 2.0) * 0.5));

  vis_nodes[0][0] += art_vis;
  vis_nodes[1][0] += art_vis;
  vis_nodes[2][0] += art_vis;
  vis_nodes[3][0] += art_vis;

  count_nodes[0][0] += 1.0;
  count_nodes[1][0] += 1.0;
  count_nodes[2][0] += 1.0;
  count_nodes[3][0] += 1.0;
}
