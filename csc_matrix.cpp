#pragma once

#include <cstdlib>
#include <cstring>

#include "matrix_utils.cpp"
#include "pcg.h"
#include "spmv_def.h"

void ldu_to_csc(const LduMatrix &ldu_matrix, CscMatrix &csc_matrix) {
    csc_matrix.cols = ldu_matrix.cells;
    csc_matrix.data_size = count_ldu_matrix_nonzero_elements(ldu_matrix);
    csc_matrix.col_off = (int *)malloc((csc_matrix.cols + 1) * sizeof(int));
    csc_matrix.rows = (int *)malloc(csc_matrix.data_size * sizeof(int));
    csc_matrix.data = (double *)malloc(csc_matrix.data_size * sizeof(double));

    int row, col, offset;
    int *tmp = (int *)malloc((csc_matrix.cols + 1) * sizeof(int));

    csc_matrix.col_off[0] = 0;
    for (int i = 1; i < ldu_matrix.cells + 1; i++) {
        if (ldu_matrix.diag[i - 1] != 0.0) {
            csc_matrix.col_off[i] = 1;
        } else {
            csc_matrix.col_off[i] = 0;
        }
    }

    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.lower[i] != 0.0) {
            col = ldu_matrix.lPtr[i];
            csc_matrix.col_off[col + 1]++;
        }
    }

    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.upper[i] != 0.0) {
            row = ldu_matrix.uPtr[i];
            csc_matrix.col_off[row + 1]++;
        }
    }

    for (int i = 0; i < ldu_matrix.cells; i++) {
        csc_matrix.col_off[i + 1] += csc_matrix.col_off[i];
    }

    memcpy(
        &tmp[0],
        &csc_matrix.col_off[0],
        (ldu_matrix.cells + 1) * sizeof(int));

    // lower
    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.lower[i] != 0.0) {
            row = ldu_matrix.uPtr[i];
            col = ldu_matrix.lPtr[i];
            offset = tmp[col]++;
            csc_matrix.rows[offset] = row;
            csc_matrix.data[offset] = ldu_matrix.lower[i];
        }
    }

    // diag
    for (int i = 0; i < ldu_matrix.cells; i++) {
        if (ldu_matrix.diag[i] != 0.0) {
            offset = tmp[i]++;
            csc_matrix.rows[offset] = i;
            csc_matrix.data[offset] = ldu_matrix.diag[i];
        }
    }

    // upper
    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.upper[i] != 0.0) {
            row = ldu_matrix.lPtr[i];
            col = ldu_matrix.uPtr[i];
            offset = tmp[col]++;
            csc_matrix.rows[offset] = row;
            csc_matrix.data[offset] = ldu_matrix.upper[i];
        }
    }

    free(tmp);
}

// basic spmv, 需要负载均衡
void csc_spmv(const CscMatrix &csc_matrix, double *vec, double *result) {
    for (int i = 0; i < csc_matrix.cols; ++i) {
        result[i] = 0;
    }

    for (int i = 0; i < csc_matrix.cols; i++) {
        int start = csc_matrix.col_off[i];
        int num = csc_matrix.col_off[i + 1] - csc_matrix.col_off[i];
        for (int j = 0; j < num; j++) {
            int row = csc_matrix.rows[start + j];
            int col = i;
            double data = csc_matrix.data[start + j];
            result[row] += vec[col] * data;
        }
    }
}

// void csc_precondition_spmv(
//     const CscMatrix &csc_matrix,
//     double *vec,
//     double *val,
//     double *result) {
//     for (int i = 0; i < csr_matrix.rows; i++) {
//         int start = csr_matrix.row_off[i];
//         int num = csr_matrix.row_off[i + 1] - csr_matrix.row_off[i];
//         double temp = 0;
//         for (int j = 0; j < num; j++) {
//             temp += vec[csr_matrix.cols[start + j]] * val[start + j];
//         }
//         result[i] = temp;
//     }
// }
//
// // diagonal precondition, get matrix M^(-1) (diagonal matrix)
// // pre_mat_val: 非对角元     : csr_matrix中元素
// //              对角元素     : 0
// // preD       : csr_matrix中对角元素的倒数
// void pcg_init_precondition_csc(const CscMatrix &csc_matrix, Precondition &pre) {
//     for (int i = 0; i < csr_matrix.rows; i++) {
//         for (int j = csr_matrix.row_off[i]; j < csr_matrix.row_off[i + 1];
//              j++) {
//             // get diagonal matrix
//             if (csr_matrix.cols[j] == i) {
//                 pre.pre_mat_val[j] = 0.;
//                 pre.preD[i] = 1.0 / csr_matrix.data[j];
//             } else {
//                 pre.pre_mat_val[j] = csr_matrix.data[j];
//             }
//         }
//     }
// }
//
// // ? 存疑，循环用处?
// void pcg_precondition_csc(
//     const CscMatrix &csc_matrix,
//     const Precondition &pre,
//     double *rAPtr,
//     double *wAPtr) {
//     double *gAPtr = (double *)malloc(csr_matrix.rows * sizeof(double));
//     v_dot_product(csr_matrix.rows, pre.preD, rAPtr, wAPtr);
//     memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
//     for (int deg = 1; deg < 2; deg++) {
//         // gAPtr = wAptr * pre.pre_mat_val; vec[rows] = matrix * vec[rows]
//         csr_precondition_spmv(csr_matrix, wAPtr, pre.pre_mat_val, gAPtr);
//         v_sub_dot_product(csr_matrix.rows, rAPtr, gAPtr, pre.preD, wAPtr);
//         memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
//     }
//     free(gAPtr);
// }

void free_csc_matrix(CscMatrix &mtx) {
    free(mtx.rows);
    free(mtx.data);
    free(mtx.col_off);
}

#define COL_SIZE(idx) (mtx.col_off[(idx) + 1] - mtx.col_off[idx])

/**
 * csc_matrix_chunk_num() - 计算给定矩阵按照每个 chunk 大小不超过 max_chunk_size 的限制条件切分时，得到的 chunk 数量。
 * 
 * @mtx: 待切分矩阵。
 * @max_chunk_size: chunk 大小的最大值限制。
 * 
 * Return: chunk 数量。
 *         若不存在一种切分方式满足每个 chunk 大小不超过 max_chunk_size 的限制，则返回 0x3f3f3f3f。
 */
int csc_matrix_chunk_num(const CscMatrix &mtx, int max_chunk_size) {
    for (int i = 0; i < mtx.cols; ++i) {
        if (COL_SIZE(i) > max_chunk_size) {
            return 0x3f3f3f3f;
        }
    }

    int chunk_cnt = 0;
    for (int i = 0, current_chunk_size = 0; i < mtx.cols; ++i) {
        current_chunk_size += COL_SIZE(i);
        if (current_chunk_size > max_chunk_size) {
            current_chunk_size = COL_SIZE(i);
            chunk_cnt += 1;
        }
    }
    chunk_cnt += 1;
    return chunk_cnt;
}

/**
 * csc_matrix_chunk_size() - 计算将给定矩阵尽可能均匀切分为若干个 chunk 时，切分后最大的 chunk 的大小。
 * 
 * @mtx: 待切分矩阵。
 * @chunk_num: 需要切分得到的 chunk 的数量。
 * 
 * Return: 切分后最大的 chunk 的大小
 */
int csc_matrix_chunk_size(const CscMatrix &mtx, int chunk_num) {
    int l = 0, r = mtx.data_size;
    while (l < r) {
        int chunk_size = l + (r - l) / 2;
        if (csc_matrix_chunk_num(mtx, chunk_size) <= chunk_num) {
            r = chunk_size;
        } else {
            l = chunk_size + 1;
        }
    }

    return l;
}

/**
 * csc_chunk_block_num() - 计算给定 chunk 中包含的 block 数量。
 * 
 * @mtx: chunk 所属矩阵。
 * @chunk_begin: chunk 的起始列号。
 * @chunk_end: chunk 的终止列号（左闭右开）。
 * 
 * 本函数会以内部的 MAX_CELL_NUM_ON_SLAVE_CORE 常量为 block 大小上限，将 chunk 切分为若干个 block。
 * 该常量的值为单个 slave 节点上能放下 double 类型数量的最大值。
 * 保证每个 block 的大小小于，并且尽量接近 MAX_CELL_NUM_ON_SLAVE_CORE。
 * 本函数会考虑 spmv 时会用到的向量切片、结果数组大小。
 * 
 * Return: block 数量。
 */
int csc_trunk_block_num(const CscMatrix &mtx, int chunk_begin, int chunk_end) {
    const int MAX_CELL_NUM_ON_SLAVE_CORE = 30000;

    int available_cell_num = MAX_CELL_NUM_ON_SLAVE_CORE;
    available_cell_num -= mtx.cols;  // 结果矩阵
    available_cell_num -= chunk_end - chunk_begin;  // 向量切片
    // TODO
}

/**
 * split_csc_matrix() - 将给定的 CscMatrix 尽可能均匀地切分成若干个 chunk
 * 
 * @mtx: 待切分矩阵
 * @slice_num: 需要切分出来的 chunk 数量。
 * 
 * Return: 数组。里面包含了 slice_num 个指针，指向切分出来的 CscChunk 结构。
 */
CscChunk **split_csc_matrix(const CscMatrix &mtx, int slice_num) {
    int max_chunk_size = csc_matrix_chunk_size(mtx, slice_num);
    int chunk_num = csc_matrix_chunk_num(mtx, max_chunk_size);

    CscChunk **result = (CscChunk **)malloc(sizeof(CscChunk *) * chunk_num);

    for (int current_chunk_size = 0,
             current_chunk_begin = 0,
             current_chunk_end = 0,
             current_chunk_idx = 0;
         current_chunk_end <= mtx.cols;
         ++current_chunk_end) {
        if (current_chunk_end != mtx.cols) {
            current_chunk_size += COL_SIZE(current_chunk_end);
        }
        if (current_chunk_size > max_chunk_size
            || current_chunk_end == mtx.cols) {
            int block_num = csc_trunk_block_num(
                mtx,
                current_chunk_begin,
                current_chunk_end);

            current_chunk_size = COL_SIZE(current_chunk_end);
            current_chunk_begin = current_chunk_end;
            current_chunk_idx += 1;
        }
    }

    return result;
}
