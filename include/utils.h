#ifndef __INS_UTILS_H
#define __INS_UTILS_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

DG_FP *getOP2PtrDevice(op_dat dat, op_access acc);
void releaseOP2PtrDevice(op_dat dat, op_access acc, const DG_FP *ptr);
DG_FP *getOP2PtrHost(op_dat dat, op_access acc);
void releaseOP2PtrHost(op_dat dat, op_access acc, const DG_FP *ptr);

DG_FP *getOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc);
void releaseOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc, const DG_FP *ptr);
DG_FP *getOP2PtrHostMap(op_dat dat, op_map map, op_access acc);
void releaseOP2PtrHostMap(op_dat dat, op_map map, op_access acc, const DG_FP *ptr);

#if defined(OP2_DG_CUDA) && DG_DIM == 3
void transfer_kernel_ptrs();
#endif

#endif
