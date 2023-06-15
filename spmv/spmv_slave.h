#ifndef _SPMV_H_
#define _SPMV_H_

#include "spmv_def.h"
#include "slave_def.h"
#include <crts.h>
#include <slave.h>
#include <string.h>

// 在从核上计算 csc 格式的稀疏矩阵乘向量(spmv)
void slave_csc_spmv(SpmvPara *para_mp);

#endif