inline void poisson_get_num_unknowns(const int *p, int *unknowns) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];
  *unknowns += dg_np;
}
