inline void ins_3d_shock_cap_art_vis_1(DG_FP *max_vis, const DG_FP *max_diff_,
                              const DG_FP *smooth_tol, const DG_FP *discon_tol,
                              const DG_FP **vis_nodes, const int **node_count,
                              DG_FP *out_vis) {
  const DG_FP max_diff = *max_diff_;
  DG_FP vertex_values[4];
  for(int i = 0; i < 4; i++) {
    vertex_values[i] = node_count[i][0] > 0 ? vis_nodes[i][0] / (DG_FP)node_count[i][0] : 0.0;
  }

  const DG_FP *r = dg_r_kernel + (DG_ORDER - 1) * DG_NP;
  const DG_FP *s = dg_s_kernel + (DG_ORDER - 1) * DG_NP;
  const DG_FP *t = dg_t_kernel + (DG_ORDER - 1) * DG_NP;

  for(int i = 0; i < DG_NP; i++) {
    out_vis[i] = 0.5 * (-(1.0 + r[i] + s[i] + t[i]) * vertex_values[0] + (1.0 + r[i]) * vertex_values[1] + (1.0 + s[i]) * vertex_values[2] + (1.0 + t[i]) * vertex_values[3]);

    if(out_vis[i] < *smooth_tol)
      out_vis[i] = 0.0;
    else if(out_vis[i] > *discon_tol)
      out_vis[i] = max_diff;
    else
      out_vis[i] = max_diff * ((out_vis[i] - *smooth_tol) / (*discon_tol - *smooth_tol));

    *max_vis = fmax(*max_vis, out_vis[i]);
  }
}