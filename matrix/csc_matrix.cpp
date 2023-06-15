#include <algorithm>
#include <cstdlib>
#include <cstring>

#include "matrix_utils.h"
#include "pcg.h"
#include "spmv_def.h"
#include "csc_matrix.h"

/**
 * csc_matrix_sort_elements() - 对 CscMatrix 中的元素排序，使每列中的元素的行号都单调递增。
 *
 * @mtx: 待排序的 CscMatrix
 */
static void csc_matrix_sort_elements(CscMatrix &mtx) {
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

extern void ldu_to_csc(const LduMatrix &ldu_matrix, CscMatrix &csc_matrix) {
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

extern void free_csc_matrix(CscMatrix &mtx) {
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
static int csc_matrix_chunk_num(const CscMatrix &mtx, int max_chunk_size) {
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
static int csc_matrix_chunk_size(const CscMatrix &mtx, int chunk_num) {
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
 * pjz 想看的东西
 */
static void split_csc_matrix_get_chunk_ranges(
    const CscMatrix &mtx,
    SplitedCscMatrix &result,
    int max_chunk_size,
    int chunk_num) {
    int i, current_chunk_size, current_chunk_begin, current_chunk_idx;
    for (i = 0,
        current_chunk_size = 0,
        current_chunk_begin = 0,
        current_chunk_idx = 0;
         i < mtx.cols;
         ++i) {
        if (current_chunk_size + COL_SIZE(i) > max_chunk_size) {
            result.chunk_ranges[current_chunk_idx].col_begin =
                current_chunk_begin;

            result.chunk_ranges[current_chunk_idx].col_num =
                i - current_chunk_begin;
            result.chunk_ranges[current_chunk_idx].size = current_chunk_size;

            current_chunk_size = 0;
            current_chunk_begin = i;
            current_chunk_idx += 1;
        }

        current_chunk_size += COL_SIZE(i);
    }

    if (current_chunk_idx < chunk_num) {
        result.chunk_ranges[current_chunk_idx].col_begin = current_chunk_begin;
        result.chunk_ranges[current_chunk_idx].col_num =
            i - current_chunk_begin;
        result.chunk_ranges[current_chunk_idx].size = current_chunk_size;

        current_chunk_size = 0;
        current_chunk_begin = i;
        current_chunk_idx += 1;
    }

    while (current_chunk_idx < chunk_num) {
        result.chunk_ranges[current_chunk_idx].col_begin = mtx.cols;
        result.chunk_ranges[current_chunk_idx].col_num = 0;
        result.chunk_ranges[current_chunk_idx].size = 0;
        current_chunk_idx += 1;
    }
}

static void build_single_chunk(
    const CscMatrix &mtx,
    CscChunk *&result,
    int result_chunk_idx,
    const CscChunkRange *ranges,
    int chunk_num) {
    int block_num = chunk_num;
    int chunk_mem_size = sizeof(CscChunk) + sizeof(CscBlock) * block_num;

    result = (CscChunk *)malloc(chunk_mem_size);
    result->size = chunk_mem_size;
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

    if (ranges[result_chunk_idx].col_num == 0) {
        return;
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
            double data = mtx.data[j];
            for (int k = 0; k < block_num; ++k) {
                if (result->blocks[k].row_begin <= row
                    && row < result->blocks[k].row_end) {
                    int &col_off =
                        result->blocks[k].col_off[i - result->col_begin];
                    result->blocks[k].data[col_off] = data;
                    result->blocks[k].rows[col_off] =
                        row - result->blocks[k].row_begin;
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

extern void csc_chunk_pack(CscChunk &chunk) {
    int meta_seg_num = 2;
    int rows_seg_num = 0;
    int col_off_seg_num = 0;
    int data_seg_num = 0;

    for (int i = 0; i < chunk.block_num; ++i) {
        CscBlock &block = chunk.blocks[i];
        data_seg_num += block.data_size;
        col_off_seg_num += block.col_num + 1;
        rows_seg_num += block.data_size;
    }

    int int_num = meta_seg_num + col_off_seg_num + rows_seg_num;
    int double_num = data_seg_num;
    chunk.packed_data.mem_size =
        sizeof(double) * (double_num + (int_num + 1) / 2);
    chunk.packed_data.mem = malloc(chunk.packed_data.mem_size);

    int *int_mem = (int *)chunk.packed_data.mem;
    double *double_mem = ((double *)chunk.packed_data.mem) + (int_num + 1) / 2;

    int *meta_seg = int_mem;
    meta_seg[0] = meta_seg_num + rows_seg_num;
    meta_seg[1] = double_mem - (double *)chunk.packed_data.mem;

    int *rows_seg = (int *)chunk.packed_data.mem + meta_seg_num;
    int *col_off_seg = (int *)chunk.packed_data.mem + meta_seg[0];
    double *data_seg = (double *)chunk.packed_data.mem + meta_seg[1];

    for (int i = 0; i < chunk.block_num; ++i) {
        CscBlock &block = chunk.blocks[i];

        memcpy(rows_seg, block.rows, sizeof(int) * block.data_size);
        rows_seg += block.data_size;
        memcpy(data_seg, block.data, sizeof(double) * block.data_size);
        data_seg += block.data_size;
        memcpy(col_off_seg, block.col_off, sizeof(int) * (block.col_num + 1));
        col_off_seg += block.col_num + 1;

        free(block.rows);
        free(block.data);
        free(block.col_off);
    }
}

extern void csc_chunk_unpack(CscChunk &chunk) {
    int meta_seg_num = 2;
    int *meta_seg = (int *)chunk.packed_data.mem;
    int *rows_seg = (int *)chunk.packed_data.mem + meta_seg_num;
    int *col_off_seg = (int *)chunk.packed_data.mem + meta_seg[0];
    double *data_seg = (double *)chunk.packed_data.mem + meta_seg[1];

    for (int i = 0; i < chunk.block_num; ++i) {
        CscBlock &block = chunk.blocks[i];

        block.rows = rows_seg;
        rows_seg += block.data_size;

        block.data = data_seg;
        data_seg += block.data_size;

        block.col_off = col_off_seg;
        col_off_seg += block.col_num + 1;
    }
}

/**
 * vaaandark 想看的东西
 */
static void split_csc_matrix_build_chunks(
    const CscMatrix &mtx,
    SplitedCscMatrix &result) {
    for (int i = 0; i < result.chunk_num; ++i) {
        build_single_chunk(
            mtx,
            result.chunks[i],
            i,
            result.chunk_ranges,
            result.chunk_num);
    }
}

extern void free_packed_splited_csc_matrix_chunk(CscChunk *chunk) {
    free(chunk->packed_data.mem);
    free(chunk);
}

extern void free_packed_splited_csc_matrix(SplitedCscMatrix *splited) {
    for (int i = 0; i < splited->chunk_num; ++i) {
        free_packed_splited_csc_matrix_chunk(splited->chunks[i]);
    }
    free(splited->chunks);
    free(splited->chunk_ranges);
}

/**
 * split_csc_matrix() - 将给定的 CscMatrix 尽可能均匀地切分成若干个 chunk
 *
 * @mtx: 待切分矩阵
 * @slice_num: 需要切分出来的 chunk 数量。
 *
 * Return: 数组。里面包含了 slice_num 个指针，指向切分出来的 CscChunk 结构。
 */
extern SplitedCscMatrix split_csc_matrix(const CscMatrix &mtx, int chunk_num) {
    SplitedCscMatrix result;
    result.chunk_num = chunk_num;
    result.chunk_ranges =
        (CscChunkRange *)malloc(sizeof(CscChunkRange) * result.chunk_num);
    result.chunks = (CscChunk **)malloc(sizeof(CscChunk *) * result.chunk_num);

    int max_chunk_size = csc_matrix_chunk_size(mtx, chunk_num);

    split_csc_matrix_get_chunk_ranges(mtx, result, max_chunk_size, chunk_num);
    split_csc_matrix_build_chunks(mtx, result);

    return result;
}
