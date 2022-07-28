inline void sample_interface(const double *x, const double *y, const double *surface,
                             double *sample_x, double *sample_y) {
  const double node0s = surface[0];
  const double node1s = surface[2];
  const double node2s = surface[5];

  double end0x = NAN;
  double end0y = NAN;
  double end1x = NAN;
  double end1y = NAN;

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

  double dist = sqrt((end1x - end0x) * (end1x - end0x) + (end1y - end0y) * (end1y - end0y));
  double dist_per_sample = dist / (LS_SAMPLE_NP + 1.0);
  
  double incrementx = ((end1x - end0x) / dist) * dist_per_sample;
  double incrementy = ((end1y - end0y) / dist) * dist_per_sample;

  for(int i = 0; i < LS_SAMPLE_NP; i++) {
    sample_x[i] = end0x + incrementx * (i + 1);
    sample_y[i] = end0y + incrementy * (i + 1);
  }
}
