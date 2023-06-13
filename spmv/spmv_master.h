#ifndef _SPMV_MASTER_H_
#define _SPMV_MASTER_H_

#include "pcg_def.h"
#include "spmv_def.h"

/// 并行计算稀疏矩阵乘向量(spmv)
void csc_spmv(const CscMatrix *csc_matrix, const double *vec, double *result);

#endif