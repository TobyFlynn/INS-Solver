inline void ins_2d_shock_cap_art_vis_1(DG_FP *max_vis, const DG_FP *max_diff_,
                              const DG_FP *smooth_tol, const DG_FP *discon_tol,
                              const DG_FP **vis_nodes, const DG_FP **node_count,
                              DG_FP *out_vis) {
  const DG_FP max_diff = *max_diff_;
  DG_FP vertex_values[3];
  for(int i = 0; i < 3; i++) {
    vertex_values[i] = node_count[i][0] > 0.0 ? vis_nodes[i][0] / node_count[i][0] : 0.0;
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

    *max_vis = fmax(*max_vis, out_vis[i]);
  }
}
