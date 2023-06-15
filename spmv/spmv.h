#ifndef _SPMV_H_
#define _SPMV_H_

/// 并行计算稀疏矩阵乘向量(spmv)
/// 包装并隐藏了在从核中计算和在主核中的规约
///
/// 使用方法：
/// ```c
/// include "spmv/spmv.h"
/// ```

#include "spmv_master.h"
#include "spmv_slave.h"

#endif