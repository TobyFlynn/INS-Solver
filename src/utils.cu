#include "op_seq.h"

#include <memory>

#include "dg_utils.h"

#define THREADS_PER_BLOCK 256

double *getOP2Array(op_dat dat) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_SUB_CELLS, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 1, args);
  double *res = (double *)malloc(dat->set->size * DG_SUB_CELLS * sizeof(double));
  cudaMemcpy(res, dat->data_d, dat->set->size * DG_SUB_CELLS * sizeof(double), cudaMemcpyDeviceToHost);
  op_mpi_set_dirtybit_cuda(1, args);
  return res;
}

__device__ double jacobiP_u(const double x, const double alpha, const double beta,
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

__device__ double gradJacobiP_u(const double x, const double alpha,
                              const double beta, const int N) {
  if(N == 0) {
    return 0.0;
  } else {
    double fact = sqrt(N * (N + alpha + beta + 1.0));
    return fact * jacobiP_u(x, alpha + 1.0, beta + 1.0, N - 1);
  }
}

__device__ double grad2JacobiP_u(const double x, const double alpha,
                               const double beta, const int N) {
  if(N == 0 || N == 1) {
    return 0.0;
  } else {
    double fact = sqrt(N * (N + alpha + beta + 1.0));
    return fact * gradJacobiP_u(x, alpha + 1.0, beta + 1.0, N - 1);
  }
}

// Calculate 2D orthonomal poly on simplex of order i,j
__device__ double simplex2DP_u(const double a, const double b, const int i,
                             const int j) {
  double h1 = jacobiP_u(a, 0, 0, i);
  double h2 = jacobiP_u(b, 2 * i + 1, 0, j);
  return sqrt(2.0) * h1 * h2 * pow(1.0 - b, i);
}

// Calculate derivatives of modal basis on simplex
__device__ void gradSimplex2DP_u(const double a, const double b, const int i,
                               const int j, double &dr, double &ds) {
  double fa  = jacobiP_u(a, 0.0, 0.0, i);
  double gb  = jacobiP_u(b, 2.0 * i + 1.0, 0.0, j);
  double dfa = gradJacobiP_u(a, 0.0, 0.0, i);
  double dgb = gradJacobiP_u(b, 2.0 * i + 1.0, 0.0, j);

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

// Calculate simplexes for Hessian
__device__ void hessianSimplex2DP_u(const double a, const double b, const int i,
                                  const int j, double &dr2, double &drs,
                                  double &ds2) {
  double fa   = jacobiP_u(a, 0.0, 0.0, i);
  double gb   = jacobiP_u(b, 2.0 * i + 1.0, 0.0, j);
  double dfa  = gradJacobiP_u(a, 0.0, 0.0, i);
  double dgb  = gradJacobiP_u(b, 2.0 * i + 1.0, 0.0, j);
  double dfa2 = grad2JacobiP_u(a, 0.0, 0.0, i);
  double dgb2 = grad2JacobiP_u(b, 2.0 * i + 1.0, 0.0, j);

  // dr2
  dr2 = dfa2 * gb;
  if(i > 1) {
    dr2 = 4.0 * dr2 * pow(1.0 - b, i - 2);
  }

  // dsr
  // dsr = dfa % gb;
  // if(i > 1) {
  //   dsr = 2.0 * dsr % arma::pow(1.0 - b, i - 2);
  // }
  //
  // arma::vec tmp = arma::vec(a.n_elem, arma::fill::zeros);
  // if(i > 0) {
  //   tmp = dgb % arma::pow(1.0 - b, i - 1);
  // }
  // if(i > 1) {
  //   tmp = tmp - i * gb % arma::pow(1.0 - b, i - 2);
  // }
  // tmp = tmp % dfa * 2.0;
  // dsr = dsr + tmp;

  drs = 2.0 * (1.0 + a) * dfa2 * gb;
  if(i > 1) {
    drs = drs * pow(1.0 - b, i - 2);
  }
  double tmp = dfa * gb * 2.0;
  if(i > 1) {
    tmp = tmp * pow(1.0 - b, i - 2);
  }
  drs = drs + tmp;
  tmp = dgb;
  if(i > 0) {
    tmp = tmp * pow(1.0 - b, i - 1);
    if(i > 1) {
      tmp = tmp - i * gb * pow(1.0 - b, i - 2);
    } else {
      tmp = tmp - i * gb;
    }
  }
  tmp = tmp * dfa * 2.0;
  drs = drs + tmp;

  // ds2
  ds2 = dfa2 * gb * pow(1.0 + a, 2);
  if(i > 1) {
    ds2 = ds2 * pow(1.0 - b, i - 2);
  }

  double tmp2 = dgb2 * pow(1.0 - b, i);
  if(i > 0) {
    tmp2 = tmp2 - 2.0 * i * dgb * pow(1.0 - b, i - 1);
  }
  if(i > 1) {
    tmp2 = tmp2 + i * (i - 1) * gb * pow(1.0 - b, i - 2);
  }
  ds2 = ds2 + fa * tmp2;

  dr2 = pow(2.0, 0.5) * dr2;
  drs = pow(2.0, 0.5) * drs;
  ds2 = pow(2.0, 0.5) * ds2;
}

__device__ double val_at_pt_u(const double r, const double s, const double *modal) {
  double a = -1.0;
  if(s != 1.0)
    a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
  double b = s;

  double new_val = 0.0;
  int modal_ind = 0;
  for(int x_ = 0; x_ <= 3; x_++) {
    for(int y_ = 0; y_ <= 3 - x_; y_++) {
      new_val += modal[modal_ind++] * simplex2DP_u(a, b, x_, y_);
    }
  }

  return new_val;
}

__device__ void grad_at_pt_u(const double r, const double s, const double *modal,
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
      gradSimplex2DP_u(a, b, x_, y_, dr_tmp, ds_tmp);
      dr += modal[modal_ind] * dr_tmp;
      ds += modal[modal_ind++] * ds_tmp;
    }
  }
}

// Get the Hessian at a point within a cell from modal values
__device__ void hessian_at_pt_u(const double r, const double s, const double *modal,
                              double &dr2, double &drs, double &ds2) {
  double a = -1.0;
  if(s != 1.0)
    a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
  double b = s;

  dr2 = 0.0;
  drs = 0.0;
  ds2 = 0.0;
  int modal_ind = 0;
  for(int x_ = 0; x_ <= 3; x_++) {
    for(int y_ = 0; y_ <= 3 - x_; y_++) {
      double dr2_tmp, drs_tmp, ds2_tmp;
      hessianSimplex2DP_u(a, b, x_, y_, dr2_tmp, drs_tmp, ds2_tmp);
      dr2 += modal[modal_ind] * dr2_tmp;
      drs += modal[modal_ind] * drs_tmp;
      ds2 += modal[modal_ind++] * ds2_tmp;
    }
  }
}

__device__ void global_xy_to_rs_u(const double x, const double y, double &r,
                                double &s, const double *cellX,
                                const double *cellY) {
  double l2 = (cellY[1] - cellY[2]) * (x - cellX[2]) + (cellX[2] - cellX[1]) * (y - cellY[2]);
  l2 = l2 / ((cellY[1] - cellY[2]) * (cellX[0] - cellX[2]) + (cellX[2] - cellX[1]) * (cellY[0] - cellY[2]));
  double l3 = (cellY[2] - cellY[0]) * (x - cellX[2]) + (cellX[0] - cellX[2]) * (y - cellY[2]);
  l3 = l3 / ((cellY[1] - cellY[2]) * (cellX[0] - cellX[2]) + (cellX[2] - cellX[1]) * (cellY[0] - cellY[2]));
  double l1 = 1.0 - l2 - l3;
  s = 2.0 * l1 - 1.0;
  r = 2.0 * l3 - 1.0;
}

__device__ void swap(double &a, double &b) {
  double c = a;
  a = b;
  b = c;
}

__device__ bool newtoncp_gepp(double A[][3], double *b) {
  for(int i = 0; i < 3; ++i) {
    int j = i;
    for(int k = i + 1; k < 3; ++k)
      if (fabs(A[k][i]) > fabs(A[j][i]))
        j = k;
    if (j != i) {
      for (int k = 0; k < 3; ++k)
        swap(A[i][k], A[j][k]);
      swap(b[i], b[j]);
    }

    if (fabs(A[i][i]) < 1.0e-8)
      return false;

    double fac = 1.0 / A[i][i];
    for (int j = i + 1; j < 3; ++j)
      A[j][i] *= fac;

    for (int j = i + 1; j < 3; ++j) {
      for (int k = i + 1; k < 3; ++k)
        A[j][k] -= A[j][i]*A[i][k];
      b[j] -= A[j][i]*b[i];
    }
  }

  for (int i = 3 - 1; i >= 0; --i) {
    double sum = 0.0;
    for (int j = i + 1; j < 3; ++j)
      sum += A[i][j]*b[j];
    b[i] = (b[i] - sum) / A[i][i];
  }

  return true;
}

__global__ void newton_kernel(const int numPts, const int *cell_ind, double *s,
                              const double *x, const double *y,
                              double *closest_x, double *closest_y,
                              const double *s_modal, const double *cell_x,
                              const double *cell_y, const double *rx,
                              const double *sx, const double *ry,
                              const double *sy, const double h) {
  // Assume only have threads in x dim
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if(tid >= numPts) return;

  const int ind_np = cell_ind[tid] * DG_NP;
  const int ind_3  = cell_ind[tid] * 3;

  const double *s_modal_l = s_modal + ind_np;
  const double *cell_x_l = cell_x + ind_3;
  const double *cell_y_l = cell_y + ind_3;

  const double x_l         = x[tid];
  const double y_l         = y[tid];
  const double closest_x_l = closest_x[tid];
  const double closest_y_l = closest_y[tid];
  const double rx_l        = rx[ind_np];
  const double sx_l        = sx[ind_np];
  const double ry_l        = ry[ind_np];
  const double sy_l        = sy[ind_np];

  double lambda = 0.0;
  double node_r, node_s, pt_r, pt_s;
  global_xy_to_rs_u(x_l, y_l, node_r, node_s, cell_x_l, cell_y_l);
  global_xy_to_rs_u(closest_x_l, closest_y_l, pt_r, pt_s, cell_x_l, cell_y_l);
  double pt_x = closest_x_l;
  double pt_y = closest_y_l;
  double init_x = closest_x_l;
  double init_y = closest_y_l;
  for(int step = 0; step < 10; step++) {
    global_xy_to_rs_u(pt_x, pt_y, pt_r, pt_s, cell_x_l, cell_y_l);
    double pt_x_old = pt_x;
    double pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    double surface = val_at_pt_u(pt_r, pt_s, s_modal_l);
    double surface_dr, surface_ds;
    grad_at_pt_u(pt_r, pt_s, s_modal_l, surface_dr, surface_ds);
    const double surface_dx = rx_l * surface_dr + sx_l * surface_ds;
    const double surface_dy = ry_l * surface_dr + sy_l * surface_ds;
    // Evaluate Hessian
    double hessian_tmp[3];
    hessian_at_pt_u(pt_r, pt_s, s_modal_l, hessian_tmp[0], hessian_tmp[1], hessian_tmp[2]);
    double hessian[3];
    hessian[0] = rx_l * rx_l * hessian_tmp[0] + sx_l * sx_l * hessian_tmp[2];
    hessian[1] = rx_l * ry_l * hessian_tmp[1] + sx_l * sy_l * hessian_tmp[1];
    hessian[2] = ry_l * ry_l * hessian_tmp[0] + sy_l * sy_l * hessian_tmp[2];

    // Check if |nabla(surface)| = 0, if so then return
    double gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy;
    if(gradsqrnorm < 1e-14)
      break;

    // Init lambda at first step
    if(step == 0)
      lambda = ((x_l - pt_x) * surface_dx + (y_l - pt_y) * surface_dy) / gradsqrnorm;

    // Gradient of functional
    double gradf[3];
    gradf[0] = pt_x - x_l + lambda * surface_dx;
    gradf[1] = pt_y - y_l + lambda * surface_dy;
    gradf[1] = surface;

    // Calculate Hessian of functional
    double hessianf[3][3];
    hessianf[0][0] = 1.0 + lambda * hessian[0];
    hessianf[0][1] = lambda * hessian[1]; hessianf[1][0] = hessianf[0][1];
    hessianf[0][2] = surface_dx; hessianf[2][0] = hessianf[0][2];

    hessianf[1][1] = 1.0 + lambda * hessian[2];
    hessianf[1][2] = surface_dy; hessianf[2][1] = hessianf[1][2];

    hessianf[2][2] = 0.0;

    if(newtoncp_gepp(hessianf, gradf)) {
      // Clamp update
      double msqr = gradf[0] * gradf[0] + gradf[1] * gradf[1];
      if(msqr > h * 0.5 * h * 0.5) {
        gradf[0] = gradf[0] * 0.5 * h / sqrt(msqr);
        gradf[1] = gradf[1] * 0.5 * h / sqrt(msqr);
        gradf[2] = gradf[2] * 0.5 * h / sqrt(msqr);
      }

      // Update guess
      pt_x -= gradf[0];
      pt_y -= gradf[1];
      lambda -= gradf[2];
    } else {
      double delta1_x = (surface / gradsqrnorm) * surface_dx;
      double delta1_y = (surface / gradsqrnorm) * surface_dy;
      lambda = ((x_l - pt_x) * surface_dx + (y_l - pt_y) * surface_dy) / gradsqrnorm;
      double delta2_x = pt_x - x_l + lambda * surface_dx;
      double delta2_y = pt_y - y_l + lambda * surface_dy;
      double msqr = delta2_x * delta2_x + delta2_y * delta2_y;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
    }

    // Gone outside the element, return
    // if((init_x - pt_x) * (init_x - pt_x) + (init_y - pt_y) * (init_y - pt_y) > h * h) {
    //   pt_x = pt_x_old;
    //   pt_y = pt_y_old;
    //   break;
    // }

    // Converged, no more steps required
    if((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) < 1e-18)
      break;
  }

  double s_ = s[tid];
  bool negative = s[tid] < 0.0;
  s[tid] = (pt_x - x_l) * (pt_x - x_l) + (pt_y - y_l) * (pt_y - y_l);
  s[tid] = sqrt(s[tid]);
  if(negative) s[tid] *= -1.0;
  if(fabs(s_ - s[tid]) < 1e-5) s[tid] = s_;
}

void newton_method(const int numPts, double *s, double *closest_x,
                   double *closest_y, int *cell_ind, const double *x,
                   const double *y, const double *s_modal, const double *cell_x,
                   const double *cell_y, const double *rx, const double *sx,
                   const double *ry, const double *sy, const double h) {
  int threadsPerBlock, blocks;
  if(numPts < THREADS_PER_BLOCK) {
    threadsPerBlock = numPts;
    blocks = 1;
  } else {
    threadsPerBlock = THREADS_PER_BLOCK;
    blocks = numPts /  THREADS_PER_BLOCK;
    if(numPts % THREADS_PER_BLOCK != 0)
      blocks++;
  }

  newton_kernel<<<blocks, threadsPerBlock>>>(numPts, cell_ind, s, x, y,
                                             closest_x, closest_y, s_modal,
                                             cell_x, cell_y, rx, sx, ry, sy, h);
}
