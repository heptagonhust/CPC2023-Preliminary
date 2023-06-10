#ifndef _SPMV_DEF_H_
#define _SPMV_DEF_H_

/// 在矩阵分块计算中，每个 block 和 vec 的一个 slice 相乘， block 的列数等于 slice 的行数
/// 同一个 chunk 中的 block 与 slice 的计算结果并起来获得 chunk_result ：
///     chunk_result = [
///                     block_result0
///                     b_1,
///                     b_2,
///                     ...,
///                     b_{block_num-1}
///                    ]
/// 最后，所有的 chunk_result 规约获得计算结果：
///     result = \sum{chunk_result}
///
/// block 中第 i 列的非零元素为 data[col_off[i]: col_off[i+1]]（左闭右开）
typedef struct {
    // block 开始和结束的行序号（左闭右开）（暂时没有用到）
    int row_begin;
    int row_end;
    int data_size;  // **这个** block 中非零元素的个数(nn0)
    double *data;   // **这个** block 中的非零元素数组，数组大小为 nn0
    int *rows;      // **这个** block 中的每个非零元素的行号，数组大小为 nn0
    int *col_off;   // **这个** block 中每列第一个非零元素在 data 数组中的索引，数组大小为 sp_col + 1
} CscBlock;

/// 在矩阵分块计算中，每个 chunk 和 vec 的一个 slice 相乘， chunk 的列数等于 slice 的行数
/// 使用 col_begin 和 col_end 来索引 slice = vec[col_begin: col_end]（左闭右开）
///
/// 分配一个 拥有 block_num 个 block 的 chunk ：
/// ```c
/// CscBlock *blk = (CscBlock *)malloc(sizeof(CscChunk) + block_num * sizeof(CscBlock))
/// ```
typedef struct {
    double *vec;   // size = sp_col
    // chunk 开始和结束的列序号（左闭右开）
    int col_begin;
    int col_end;
    // chunk 开始和结束的行序号（左闭右开）（暂时没有用到）
    int row_begin;
    int row_end;
    // 将一个 chunk 分成 block_num 个 block
    int block_num;
    CscBlock blocks[];
} CscChunk;

/// 从核计算 spmv 时传递的参数，逻辑上包含了（结构上不直接包含）：
///     - 稀疏矩阵
///     - 右乘向量
///     - 结果向量
/// 矩阵被分块到了不同的 chunk 和 block 中，右乘向量也在 chunk 中
typedef struct {
    // 将一个 sp 分成 64 个 chunk
    int chunk_num;      // 64
    CscChunk *chunks;
    int sp_row;         // 稀疏矩阵的行数
    int sp_col;         // 稀疏矩阵的列数，也是右乘向量的行数
    double *result;     // 结果列向量，行数为 sp_row
} SpmvPara;

#endif