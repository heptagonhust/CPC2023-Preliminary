#ifndef _PCG_DEF_H_
#define _PCG_DEF_H_

#include "spmv_slave_def.h"

typedef struct {
    int col_begin;
    int col_num;
    int size;
} CscChunkRange;

typedef struct {
    int chunk_num;
    CscChunkRange *chunk_ranges;
    CscChunk **chunks;
} SplitedCscMatrix;

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

/**
 * struct CsrMatrix - 压缩稀疏行矩阵
 * @rows: 行数
 * @row_off: 行偏移数组，大小为 cols + 1，每一行的第一个非零元素在 data 数组中的位置。最后一个元素等于 data_size。
 * @cols: 列索引数组，大小为 data_size，每个非元素所在的列
 * @data: 大小为 data_size，每个非零元素的值
 * @data_size: 非零元素个数
 */
struct CsrMatrix {
    int rows;
    int *row_off;
    int *cols;
    double *data;
    int data_size;
};

struct LduMatrix {
    double *upper;
    double *lower;
    double *diag;
    int *uPtr;
    int *lPtr;
    int faces;
    int cells;
};

struct PCG {
    double *r;
    double *z;
    double *p;
    double *Ax;
    double sumprod;
    double sumprod_old;
    double residual;
    double alpha;
    double beta;

    double *x;
    double *source;
};

struct Precondition {
    double *pre_mat_val;
    double *preD;
};

struct PCGReturn {
    double residual;
    int iter;
};

#include <stdio.h>
#include <time.h>
#define INFO(M, ...) \
    { \
        time_t t; \
        struct tm *tmif; \
        t = time(NULL); \
        tmif = localtime(&t); \
        printf( \
            "[%d-%2d-%2d] [%2d:%2d:%2d] [INFO] " M "", \
            tmif->tm_year + 1900, \
            tmif->tm_mon + 1, \
            tmif->tm_mday, \
            tmif->tm_hour, \
            tmif->tm_min, \
            tmif->tm_sec, \
            ##__VA_ARGS__); \
    }

#endif
