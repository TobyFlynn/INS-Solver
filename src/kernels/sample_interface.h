bool is_point_in_cell_(const double x, const double y, const double *cellX, const double *cellY) {
  double ABx = cellX[1] - cellX[0];
  double ABy = cellY[1] - cellY[0];

  double APx = x - cellX[0];
  double APy = y - cellY[0];

  double BCx = cellX[2] - cellX[1];
  double BCy = cellY[2] - cellY[1];

  double BPx = x - cellX[1];
  double BPy = y - cellY[1];

  double CAx = cellX[0] - cellX[2];
  double CAy = cellY[0] - cellY[2];

  double CPx = x - cellX[2];
  double CPy = y - cellY[2];

  double AB_AP = ABx * APy - ABy * APx;
  double BC_BP = BCx * BPy - BCy * BPx;
  double CA_CP = CAx * CPy - CAy * CPx;

  bool zero0 = AB_AP == 0.0 || AB_AP == -0.0;
  bool zero1 = BC_BP == 0.0 || BC_BP == -0.0;
  bool zero2 = CA_CP == 0.0 || CA_CP == -0.0;

  if(zero0 || zero1 || zero2) {
    // One zero means on an edge of the triangle
    // Two zeros means on a vertex of the triangle
    if(zero0 && !zero1 && !zero2) {
      return (BC_BP > 0.0 && CA_CP > 0.0) || (BC_BP < 0.0 && CA_CP < 0.0);
    } else if(!zero0 && zero1 && !zero2) {
      return (AB_AP > 0.0 && CA_CP > 0.0) || (AB_AP < 0.0 && CA_CP < 0.0);
    } else if(!zero0 && !zero1 && zero2) {
      return (AB_AP > 0.0 && BC_BP > 0.0) || (AB_AP < 0.0 && BC_BP < 0.0);
    } else if(zero0 && zero1 && !zero2) {
      return true;
    } else if(zero0 && !zero1 && zero2) {
      return true;
    } else if(!zero0 && zero1 && zero2) {
      return true;
    } else {
      return false;
    }
  }

  if(AB_AP > 0.0 && BC_BP > 0.0 && CA_CP > 0.0) {
    return true;
  } else if(AB_AP < 0.0 && BC_BP < 0.0 && CA_CP < 0.0) {
    return true;
  } else {
    return false;
  }
}

void rs_to_global_xy(const double r, const double s, double &x, double &y,
                     const double *cellX, const double *cellY) {
  x = 0.5 * (-(r + s) * cellX[0] + (1.0 + r) * cellX[1] + (1.0 + s) * cellX[2]);
  y = 0.5 * (-(r + s) * cellY[0] + (1.0 + r) * cellY[1] + (1.0 + s) * cellY[2]);
}

void global_xy_to_rs(const double x, const double y, double &r, double &s,
                     const double *cellX, const double *cellY) {
  double l2 = (cellY[1] - cellY[2]) * (x - cellX[2]) + (cellX[2] - cellX[1]) * (y - cellY[2]);
  l2 = l2 / ((cellY[1] - cellY[2]) * (cellX[0] - cellX[2]) + (cellX[2] - cellX[1]) * (cellY[0] - cellY[2]));
  double l3 = (cellY[2] - cellY[0]) * (x - cellX[2]) + (cellX[0] - cellX[2]) * (y - cellY[2]);
  l3 = l3 / ((cellY[1] - cellY[2]) * (cellX[0] - cellX[2]) + (cellX[2] - cellX[1]) * (cellY[0] - cellY[2]));
  double l1 = 1.0 - l2 - l3;
  s = 2.0 * l1 - 1.0;
  r = 2.0 * l3 - 1.0;
}

double jacobiP(const double x, const double alpha, const double beta,
               const int N) {
  double gamma0 = pow(2.0, alpha + beta + 1.0) / (alpha + beta + 1.0) *
                  tgamma(alpha + 1.0) * tgamma(beta + 1.0) /
                  tgamma(alpha + beta + 1.0);
  double p_0 = 1.0 / sqrt(gamma0);

  // First base case
  if(N == 0) {
    return p_0;
  }

  double gamma1 = (alpha + 1.0) * (beta + 1.0) / (alpha + beta + 3.0) * gamma0;
  double p_1 = ((alpha + beta + 2.0) * x / 2.0 + (alpha - beta) / 2.0) /
               sqrt(gamma1);
  // Second base case
  if(N == 1) {
    return p_1;
  }

  double p_2 = 0.0;

  // Recurrence for N > 1
  double aOld = 2.0 / (2.0 + alpha + beta) * sqrt((alpha + 1.0) * (beta + 1.0) /
                (alpha + beta + 3.0));
  for(int i = 1; i < N; i++) {
    double h1 = 2.0 * i + alpha + beta;
    double aNew = 2.0 / (h1 + 2.0) * sqrt((i + 1.0) * (i + 1.0 + alpha + beta) *
                  (i + 1.0 + alpha) * (i + 1.0 + beta) / (h1 + 1.0) /
                  (h1 + 3.0));
    double bNew = -(alpha * alpha - beta * beta) / h1 / (h1 + 2.0);
    p_2 = 1.0 / aNew * (-aOld * p_0 + (x - bNew) * p_1);
    aOld = aNew;
    p_0 = p_1;
    p_1 = p_2;
  }

  return p_2;
}

double simplex2DP(const double a, const double b, const int i, const int j) {
  double h1 = jacobiP(a, 0, 0, i);
  double h2 = jacobiP(b, 2 * i + 1, 0, j);
  return sqrt(2.0) * h1 * h2 * pow(1.0 - b, i);
}

double val_at_pt(const double r, const double s, const double *modal) {
  double a = -1.0;
  if(s != 1.0)
    a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
  double b = s;

  double new_val = 0.0;
  int modal_ind = 0;
  for(int x_ = 0; x_ <= 3; x_++) {
    for(int y_ = 0; y_ <= 3 - x_; y_++) {
      new_val += modal[modal_ind++] * simplex2DP(a, b, x_, y_);
    }
  }

  return new_val;
}

double gradJacobiP(const double x, const double alpha, const double beta,
                   const int N) {
  if(N == 0) {
    return 0.0;
  } else {
    double fact = sqrt(N * (N + alpha + beta + 1.0));
    return fact * jacobiP(x, alpha + 1.0, beta + 1.0, N - 1);
  }
}

void gradSimplex2DP(const double a, const double b, const int i, const int j,
                    double &dr, double &ds) {
  double fa  = jacobiP(a, 0.0, 0.0, i);
  double gb  = jacobiP(b, 2.0 * i + 1.0, 0.0, j);
  double dfa = gradJacobiP(a, 0.0, 0.0, i);
  double dgb = gradJacobiP(b, 2.0 * i + 1.0, 0.0, j);

  // r derivative
  dr = dfa * gb;
  if(i > 0) {
    dr = dr * pow(0.5 * (1.0 - b), i - 1);
  }

  // s derivative
  ds = dfa * (gb * (0.5 * (1.0 + a)));
  if(i > 0) {
    ds = ds * pow(0.5 * (1.0 - b), i - 1);
  }

  double tmp = dgb * pow(0.5 * (1.0 - b), i);
  if(i > 0) {
    tmp = tmp - 0.5 * i * gb * pow(0.5 * (1.0 - b), i - 1);
  }
  ds = ds + fa * tmp;

  // Normalise
  dr = pow(2.0, i + 0.5) * dr;
  ds = pow(2.0, i + 0.5) * ds;
}

void grad_at_pt(const double r, const double s, const double *modal,
                double &dr, double &ds) {
  double a = -1.0;
  if(s != 1.0)
    a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
  double b = s;

  dr = 0.0;
  ds = 0.0;
  int modal_ind = 0;
  for(int x_ = 0; x_ <= 3; x_++) {
    for(int y_ = 0; y_ <= 3 - x_; y_++) {
      double dr_tmp, ds_tmp;
      gradSimplex2DP(a, b, x_, y_, dr_tmp, ds_tmp);
      dr += modal[modal_ind] * dr_tmp;
      ds += modal[modal_ind++] * ds_tmp;
    }
  }
}

inline void sample_interface(const double *r, const double *s, const double *x,
                             const double *y, const double *surface,
                             const double *s_modal, const double *rx,
                             const double *sx, const double *ry,
                             const double *sy, double *sample_x,
                             double *sample_y) {
  // Check that the cell contains the interface
  // bool positive0 = surface[0] > 0.0;
  // bool interface = false;
  // for(int i = 1; i < DG_NP; i++) {
  //   if(positive0 != surface[i] > 0.0)
  //     interface = true;
  // }
  // if(!interface) {
  //   for(int i = 0; i < LS_SAMPLE_NP; i++) {
  //     sample_x[i] = NAN;
  //     sample_y[i] = NAN;
  //   }
  //   return;
  // }

  double sample_r[LS_SAMPLE_NP], sample_s[LS_SAMPLE_NP];

  // Initial positions of sample points (in r-s coords)
  for(int i = 0; i < LS_SAMPLE_NP; i++) {
    sample_r[i] = r[i] * 0.9;
    sample_s[i] = s[i] * 0.9;
    rs_to_global_xy(sample_r[i], sample_s[i], sample_x[i], sample_y[i], x, y);
  }

  for(int p = 0; p < LS_SAMPLE_NP; p++) {
    bool converged = false;
    for(int step = 0; step < 10; step++) {
      double surf = val_at_pt(sample_r[p], sample_s[p], s_modal);
      double dsdr, dsds;
      grad_at_pt(sample_r[p], sample_s[p], s_modal, dsdr, dsds);
      double dsdx = rx[0] * dsdr + sx[0] * dsds;
      double dsdy = ry[0] * dsdr + sy[0] * dsds;

      double sqrnorm = dsdx * dsdx + dsdy * dsdy;
      if(sqrnorm > 0.0) {
        dsdx *= surf / sqrnorm;
        dsdy *= surf / sqrnorm;
      }

      sample_x[p] -= dsdx;
      sample_y[p] -= dsdy;

      // Converged
      if(dsdx * dsdx + dsdy * dsdy < 1e-12) {
        converged = true;
        break;
      }

      global_xy_to_rs(sample_x[p], sample_y[p], sample_r[p], sample_s[p], x, y);
    }

    // Convert to x-y coords
    if(converged) {
      // Check within original element
      if(!is_point_in_cell_(sample_x[p], sample_y[p], x, y)) {
        sample_x[p] = NAN;
        sample_y[p] = NAN;
      }
    } else {
      // TODO handle this properly
      sample_x[p] = NAN;
      sample_y[p] = NAN;
    }
  }
}
