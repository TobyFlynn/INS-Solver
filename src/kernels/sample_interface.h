inline void sample_interface(const int *p, const DG_FP *x, const DG_FP *y,
                             const DG_FP *surface, DG_FP *sample_x, DG_FP *sample_y) {
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask = &FMASK[(*p - 1) * 3 * DG_NPF];

  int fmask_ind_tmp = fmask[0];
  const DG_FP node0s = surface[fmask_ind_tmp];
  fmask_ind_tmp = fmask[dg_npf - 1];
  const DG_FP node1s = surface[fmask_ind_tmp];
  fmask_ind_tmp = fmask[2 * dg_npf - 1];
  const DG_FP node2s = surface[fmask_ind_tmp];

  DG_FP end0x = NAN;
  DG_FP end0y = NAN;
  DG_FP end1x = NAN;
  DG_FP end1y = NAN;

  if((node0s > 0.0) != (node1s > 0.0)) {
    end0x = x[0] - (node0s / (node0s - node1s)) * (x[0] - x[1]);
    end0y = y[0] - (node0s / (node0s - node1s)) * (y[0] - y[1]);
  }
  if((node1s > 0.0) != (node2s > 0.0)) {
    if(isnan(end0x)) {
      end0x = x[1] - (node1s / (node1s - node2s)) * (x[1] - x[2]);
      end0y = y[1] - (node1s / (node1s - node2s)) * (y[1] - y[2]);
    } else {
      end1x = x[1] - (node1s / (node1s - node2s)) * (x[1] - x[2]);
      end1y = y[1] - (node1s / (node1s - node2s)) * (y[1] - y[2]);
    }
  }
  if((node2s > 0.0) != (node0s > 0.0)) {
    if(isnan(end0x)) {
      end0x = x[2] - (node2s / (node2s - node0s)) * (x[2] - x[0]);
      end0y = y[2] - (node2s / (node2s - node0s)) * (y[2] - y[0]);
    } else {
      end1x = x[2] - (node2s / (node2s - node0s)) * (x[2] - x[0]);
      end1y = y[2] - (node2s / (node2s - node0s)) * (y[2] - y[0]);
    }
  }

  if(isnan(end0x) || isnan(end0y) || isnan(end1x) || isnan(end1y)) {
    for(int i = 0; i < LS_SAMPLE_NP; i++) {
      sample_x[i] = NAN;
      sample_y[i] = NAN;
    }
    return;
  }

  DG_FP dist = sqrt((end1x - end0x) * (end1x - end0x) + (end1y - end0y) * (end1y - end0y));
  DG_FP dist_per_sample = dist / (LS_SAMPLE_NP + 1.0);

  DG_FP incrementx = ((end1x - end0x) / dist) * dist_per_sample;
  DG_FP incrementy = ((end1y - end0y) / dist) * dist_per_sample;

  for(int i = 0; i < LS_SAMPLE_NP; i++) {
    sample_x[i] = end0x + incrementx * (i + 1);
    sample_y[i] = end0y + incrementy * (i + 1);
  }
}
