inline void ins_2d_st_art_vis_1(const DG_FP *max_diff_, const DG_FP *smooth_tol,
                                const DG_FP *discon_tol, const DG_FP *diff_width_,
                                const DG_FP *trans_width_, const DG_FP *alpha_,
                                const DG_FP *s, const DG_FP **vis_nodes,
                                const int **node_count, DG_FP *out_vis) {
  const DG_FP PI = 3.141592653589793238463;
  const DG_FP alpha = *alpha_ * *diff_width_;
  const DG_FP max_diff = *max_diff_;
  // const DG_FP max_diff_interface = max_diff * 0.1;
  const DG_FP max_diff_interface = max_diff;
  const DG_FP trans_width = *alpha_ * *trans_width_;
  DG_FP vertex_values[3];
  for(int i = 0; i < 3; i++) {
    vertex_values[i] = node_count[i][0] > 0 ? vis_nodes[i][0] / (DG_FP)node_count[i][0] : 0.0;
  }

  const DG_FP *r_ = &dg_r_kernel[(DG_ORDER - 1) * DG_NP];
  const DG_FP *s_ = &dg_s_kernel[(DG_ORDER - 1) * DG_NP];

  for(int i = 0; i < DG_NP; i++) {
    out_vis[i]  = 0.5 * vertex_values[1] * (1.0 + r_[i]);
    out_vis[i] += 0.5 * vertex_values[2] * (1.0 + s_[i]);
    out_vis[i] -= 0.5 * vertex_values[0] * (s_[i] + r_[i]);

    if(out_vis[i] < *smooth_tol)
      out_vis[i] = 0.0;
    else if(out_vis[i] > *discon_tol)
      out_vis[i] = max_diff;
    else
      out_vis[i] = max_diff * ((out_vis[i] - *smooth_tol) / (*discon_tol - *smooth_tol));

    DG_FP tmp = fmax(1e-6, fmin(max_diff_interface, max_diff_interface * (trans_width - (fabs(s[i]) - alpha)) / trans_width));
    out_vis[i] = fmax(tmp, out_vis[i]);
/*
    if(out_vis[i] < *smooth_tol)
      out_vis[i] = 0.0;
    else if(out_vis[i] > *discon_tol)
      out_vis[i] = max_diff;
    else
      out_vis[i] = max_diff * 0.5 * (1.0 + sin(PI * (out_vis[i] - 0.5 * (*discon_tol + *smooth_tol)) / (*discon_tol - *smooth_tol)));
*/
  }
}
