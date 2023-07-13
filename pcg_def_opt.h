#pragma once

#include "pcg_def.h"
#include "spmv_def.h"

/**
 * struct CscMatrix - 压缩稀疏列矩阵
 * @cols: 列数
 * @col_off: 列偏移数组，大小为 cols + 1
 * @rows: 行索引数组，大小为 data_size
 * @data: 非零元素数组，大小为 data_size
 * @data_size: 非零元素的个数
 */
typedef struct {
    int cols;
    int data_size;
    int *col_off;
    int *rows;
    double *data;
} CscMatrix;

typedef struct {
    int row_begin;
    int row_num;
    int size;
} CooChunkRange;

typedef struct {
    int chunk_num;
    CooChunkRange *chunk_ranges;
    SizedCooChunk *chunks;
} SplitedCooMatrix;