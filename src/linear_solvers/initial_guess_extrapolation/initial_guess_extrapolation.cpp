#include "linear_solvers/initial_guess_extrapolation.h"

#include "op_seq.h"

#include "dg_utils.h"
#include "utils.h"

using namespace std;

#define INTERP_ORDER 3

void initial_guess_extrapolation(DGMesh *mesh, vector<pair<DG_FP,DGTempDat>> &history,
                                 op_dat init_guess, const DG_FP t_n1) {
  const int history_size = history.size();
  if(history.size() != EXTRAPOLATE_HISTORY_SIZE) return;

  double min_t = history[0].first;
  double max_t = history[0].first;
  for(int i = 1; i < history_size; i++) {
    min_t = history[i].first < min_t ? history[i].first : min_t;
    max_t = history[i].first > max_t ? history[i].first : max_t;
  }

  // Normalise times to [-1,1]
  const DG_FP diff_t = max_t - min_t;
  arma::vec normalised_times(history_size);
  for(int i = 0; i < history_size; i++) {
    DG_FP t = history[i].first;
    t = 2.0 * ((t - min_t) / diff_t) - 1.0;
    normalised_times[i] = t;
  }
  // std::cout << normalised_times << std::endl;

  arma::mat V_t = DGUtils::vandermonde1D(normalised_times, INTERP_ORDER);
  arma::vec one(1);
  one[0] = 2.0 * ((t_n1 - min_t) / diff_t) - 1.0;
  arma::mat V_tn1 = DGUtils::vandermonde1D(one, INTERP_ORDER);

  // arma::mat beta = V_tn1 * arma::inv(V_t.t() * V_t) * V_t.t();
  arma::vec beta = arma::solve(V_t.t(), V_tn1.t());
  // std::cout << beta << std::endl;

  DG_FP beta_arr[EXTRAPOLATE_HISTORY_SIZE];
  DG_FP beta_total = 0.0;
  for(int i = 0; i < EXTRAPOLATE_HISTORY_SIZE; i++) {
    beta_arr[i] = beta[i];
    beta_total += beta[i];
  }
  // std::cout << beta_total << std::endl;

  op_par_loop(initial_guess_extrap, "initial_guess_extrap", mesh->cells,
              op_arg_gbl(beta_arr, 8, DG_FP_STR, OP_READ),
              op_arg_dat(history[0].second.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(history[1].second.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(history[2].second.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(history[3].second.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(history[4].second.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(history[5].second.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(history[6].second.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(history[7].second.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(init_guess, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}
