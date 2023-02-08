#ifndef __INS_UTILS_H
#define __INS_UTILS_H

#include "op_seq.h"

double *getOP2PtrDevice(op_dat dat, op_access acc);
void releaseOP2PtrDevice(op_dat dat, op_access acc, const double *ptr);
double *getOP2PtrHost(op_dat dat, op_access acc);
void releaseOP2PtrHost(op_dat dat, op_access acc, const double *ptr);

double *getOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc);
void releaseOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc, const double *ptr);
double *getOP2PtrHostMap(op_dat dat, op_map map, op_access acc);
void releaseOP2PtrHostMap(op_dat dat, op_map map, op_access acc, const double *ptr);

#endif
