#ifndef _SPMV_H_
#define _SPMV_H_

#include "spmv_def.h"
#include "slave_def.h"
#include <crts.h>
#include <slave.h>
#include <string.h>

/// 在从核上计算 csc 格式的稀疏矩阵乘向量(spmv)
///
/// * `para_mp` - 从核计算所需要参数在主存中的地址

#ifdef __cplusplus
extern "C" {
#endif

  void slave_csc_spmv(SpmvPara *para_mp);

#ifdef __cplusplus
}
#endif

#endif