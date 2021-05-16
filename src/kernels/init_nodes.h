inline void init_nodes(const double **nc, double *nodeX, double *nodeY) {
  nodeX[0] = nc[0][0];
  nodeX[1] = nc[1][0];
  nodeX[2] = nc[2][0];
  nodeY[0] = nc[0][1];
  nodeY[1] = nc[1][1];
  nodeY[2] = nc[2][1];
}
