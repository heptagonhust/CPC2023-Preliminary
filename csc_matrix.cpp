#pragma once

#include <algorithm>
#include <cstdlib>
#include <cstring>

#include "matrix_utils.cpp"
#include "pcg.h"
#include "spmv_slave_def.h"

/**
 * csc_matrix_sort_elements() - 对 CscMatrix 中的元素排序，使每列中的元素的行号都单调递增。
 * 
 * @mtx: 待排序的 CscMatrix
 */
void csc_matrix_sort_elements(CscMatrix &mtx) {
    for (int i = 0; i < mtx.cols; ++i) {
        int size = mtx.col_off[i + 1] - mtx.col_off[i];
        int *rows = mtx.rows + mtx.col_off[i];
        double *data = mtx.data + mtx.col_off[i];
        struct Element {
            int row;
            double data;

            bool operator<(const Element &another) {
                return row < another.row;
            }
        };
        Element *elements = (Element *)malloc(sizeof(Element) * size);
        for (int j = 0; j < size; ++j) {
            elements[j].row = rows[j];
            elements[j].data = data[j];
        }
        std::sort(elements, elements + size);
        free(elements);
    }
}

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

    csc_matrix_sort_elements(csc_matrix);
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
 * pjz 想看的东西
 */
void split_csc_matrix_get_chunk_ranges(
    const CscMatrix &mtx,
    SplitedCscMatrix &result,
    int max_chunk_size) {
    int i, current_chunk_size, current_chunk_begin, current_chunk_idx;
    for (i = 0,
        current_chunk_size = 0,
        current_chunk_begin = 0,
        current_chunk_idx = 0;
         i < mtx.cols;
         ++i) {
        if (current_chunk_size + COL_SIZE(i) > max_chunk_size
            || i == mtx.cols - 1) {
            result.chunk_ranges[i].col_begin = current_chunk_begin;
            result.chunk_ranges[i].col_num =
                current_chunk_idx - current_chunk_begin;
            result.chunk_ranges[i].size = current_chunk_size;

            current_chunk_size = 0;
            current_chunk_begin = i;
            current_chunk_idx += 1;
        }

        current_chunk_size += COL_SIZE(i);
    }
}

void build_single_chunk(
    const CscMatrix &mtx,
    CscChunk *&result,
    int result_chunk_idx,
    const CscChunkRange *ranges,
    int chunk_num,
    double *vec) {
    int block_num = chunk_num;
    int chunk_mem_size = sizeof(CscChunk) + sizeof(CscBlock) * block_num;

    result = (CscChunk *)malloc(chunk_mem_size);
    result->size = chunk_mem_size;
    result->vec = vec + ranges[result_chunk_idx].col_begin;
    result->col_begin = ranges[result_chunk_idx].col_begin;
    result->col_end =
        ranges[result_chunk_idx].col_begin + ranges[result_chunk_idx].col_num;
    result->row_begin = 0;
    result->row_end = mtx.cols;
    result->block_num = block_num;

    for (int i = 0; i < block_num; ++i) {
        result->blocks[i].col_num = ranges[result_chunk_idx].col_num;
        result->blocks[i].row_begin = ranges[i].col_begin;
        result->blocks[i].row_end = ranges[i].col_begin + ranges[i].col_num;
        result->blocks[i].data_size = 0;
        result->blocks[i].col_off =
            (int *)malloc(sizeof(int) * (result->blocks[i].col_num + 1));
        for (int j = 0; j <= result->blocks[i].col_num; ++j) {
            result->blocks[i].col_off[j] = 0;
        }
    }

    for (int i = result->col_begin; i < result->col_end; ++i) {
        for (int j = mtx.col_off[i]; j < mtx.col_off[i + 1]; ++j) {
            int row = mtx.rows[j];
            for (int k = 0; k < block_num; ++k) {
                if (result->blocks[k].row_begin <= row
                    && row < result->blocks[k].row_end) {
                    result->blocks[k].data_size += 1;
                    result->blocks[k].col_off[i - result->col_begin] += 1;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < block_num; ++i) {
        int *col_off = result->blocks[i].col_off;
        int col_num = result->blocks[i].col_num;
        col_off[col_num] = col_off[col_num - 1];
        for (int j = 1; j <= col_num; ++j) {
            col_off[j] += col_off[j - 1];
        }
        for (int j = col_num - 1; j >= 1; --j) {
            col_off[j] -= col_off[j] - col_off[j - 1];
        }
        col_off[0] = 0;
    }

    for (int i = 0; i < block_num; ++i) {
        result->blocks[i].data =
            (double *)malloc(sizeof(double) * result->blocks[i].data_size);
        result->blocks[i].rows =
            (int *)malloc(sizeof(int) * result->blocks[i].data_size);
    }

    for (int i = result->col_begin; i < result->col_end; ++i) {
        for (int j = mtx.col_off[i]; j < mtx.col_off[i + 1]; ++j) {
            int row = mtx.rows[j];
            double data = mtx.rows[j];
            for (int k = 0; k < block_num; ++k) {
                if (result->blocks[k].row_begin <= row
                    && row < result->blocks[k].row_end) {
                    int &col_off =
                        result->blocks[k].col_off[i - result->col_begin];
                    result->blocks[k].data[col_off] = data;
                    result->blocks[k].rows[col_off] = row;
                    col_off += 1;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < block_num; ++i) {
        int *col_off = result->blocks[i].col_off;
        int col_num = result->blocks[i].col_num;
        for (int j = col_num; j >= 1; --j) {
            col_off[j] -= col_off[j] - col_off[j - 1];
        }
        col_off[0] = 0;
    }
}

/**
 * vaaandark 想看的东西
 */
void split_csc_matrix_build_chunks(
    const CscMatrix &mtx,
    SplitedCscMatrix &result,
    double *vec) {
    for (int i = 0; i < result.chunk_num; ++i) {
        build_single_chunk(
            mtx,
            result.chunks[i],
            i,
            result.chunk_ranges,
            result.chunk_num,
            vec);
    }
}

/**
 * split_csc_matrix() - 将给定的 CscMatrix 尽可能均匀地切分成若干个 chunk
 * 
 * @mtx: 待切分矩阵
 * @slice_num: 需要切分出来的 chunk 数量。
 * @vec: 将要被乘的向量
 * 
 * Return: 数组。里面包含了 slice_num 个指针，指向切分出来的 CscChunk 结构。
 */
SplitedCscMatrix
split_csc_matrix(const CscMatrix &mtx, int chunk_num, double *vec) {
    SplitedCscMatrix result;
    result.chunk_num = chunk_num;
    result.chunk_ranges =
        (CscChunkRange *)malloc(sizeof(CscChunkRange) * result.chunk_num);
    result.chunks = (CscChunk **)malloc(sizeof(CscChunk *) * result.chunk_num);

    int max_chunk_size = csc_matrix_chunk_size(mtx, chunk_num);

    split_csc_matrix_get_chunk_ranges(mtx, result, max_chunk_size);
    split_csc_matrix_build_chunks(mtx, result, vec);

    return result;
}
