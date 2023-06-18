#ifndef _SPMV_MASTER_H_
#define _SPMV_MASTER_H_

#include <stdlib.h>
#include <string.h>

#include "pcg_def.h"
#include "spmv_def.h"
#include "spmv_slave.h"
#include "csc_matrix.h"
/// 并行计算稀疏矩阵乘向量(spmv)
///
/// 将计算任务分配到所有(64个)从核进行计算
/// 
/// * `para` - 从核计算所需参数，包含了从核计算 spmv 的各种参数，用于调用 `slave_csc_spmv`
/// * `result` - spmv 的计算结果，应当是一个在调用本函数之前就分配好的数组
void coo_spmv(SpmvPara *para, double *result);

/// 将划分好的 csc 矩阵转换成 spmv_para
///
/// 需要填入的参数包括右乘的向量
///
/// * `mat` - 划分好后的 CSC 格式稀疏矩阵
/// * `para` - 待填充的从核计算所需参数
/// * `vec` - spmv 中的右乘向量
/// * `row_num` - 系数矩阵的行数
/// * `col_num` - 系数矩阵的列数
void spmv_para_from_splited_coo_matrix(
    const SplitedCooMatrix *mat,
    SpmvPara *para,
    double *vec,
    int row_num,
    int col_num);

/// 释放 spmv_para 的空间
///
/// 只释放填充 spmv_para 时分配的空间，不释放 splited_csc_matrix 分配的空间
///
/// * `para` - 需要释放的从核计算所需参数
void spmv_para_free(SpmvPara *para);

#endif