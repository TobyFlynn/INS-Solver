inline void modal_shock_detector_2d_1(const DG_FP **vis_nodes,
                                      const DG_FP **count_nodes, DG_FP *out) {
  DG_FP vertex_values[3];
  for(int i = 0; i < 3; i++) {
    vertex_values[i] = count_nodes[i][0] > 0.0 ? vis_nodes[i][0] / count_nodes[i][0] : 0.0;
  }

  const DG_FP *r_ = &dg_r_kernel[(DG_ORDER - 1) * DG_NP];
  const DG_FP *s_ = &dg_s_kernel[(DG_ORDER - 1) * DG_NP];

  for(int i = 0; i < DG_NP; i++) {
    out[i]  = 0.5 * vertex_values[1] * (1.0 + r_[i]);
    out[i] += 0.5 * vertex_values[2] * (1.0 + s_[i]);
    out[i] -= 0.5 * vertex_values[0] * (s_[i] + r_[i]);
  }
}
