#include "coo_matrix.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "matrix_utils.h"
#include "pcg.h"
#include "spmv_def.h"

/**
 * coo_matrix_chunk_num() - 计算给定矩阵按照每个 chunk 大小不超过 max_chunk_size 的限制条件切分时，得到的 chunk 数量。
 *
 * @mtx: 待切分矩阵。
 * @max_chunk_size: chunk 大小的最大值限制。
 *
 * Return: chunk 数量。
 *         若不存在一种切分方式满足每个 chunk 大小不超过 max_chunk_size 的限制，则返回 0x3f3f3f3f。
 */
static int splited_chunk_num(int *col_size, int col_num, int max_chunk_size) {
    for (int i = 0; i < col_num; ++i) {
        if (col_size[i] > max_chunk_size) {
            return 0x3f3f3f3f;
        }
    }

    int chunk_cnt = 0;
    for (int i = 0, current_chunk_size = 0; i < col_num; ++i) {
        current_chunk_size += col_size[i];
        if (current_chunk_size > max_chunk_size) {
            current_chunk_size = col_size[i];
            chunk_cnt += 1;
        }
    }

    chunk_cnt += 1;
    return chunk_cnt;
}

/**
 * coo_matrix_chunk_size() - 计算将给定矩阵尽可能均匀切分为若干个 chunk 时，切分后最大的 chunk 的大小。
 *
 * @mtx: 待切分矩阵。
 * @chunk_num: 需要切分得到的 chunk 的数量。
 *
 * Return: 切分后最大的 chunk 的大小
 */
static int splited_max_chunk_size(int *col_size, int col_num, int chunk_num) {
    int l = 0, r = 0;
    for (int i = 0; i < col_num; ++i) {
        r += col_size[i];
    }
    while (l < r) {
        int chunk_size = l + (r - l) / 2;
        if (splited_chunk_num(col_size, col_num, chunk_size) <= chunk_num) {
            r = chunk_size;
        } else {
            l = chunk_size + 1;
        }
    }

    return l;
}

static void get_chunk_ranges(
    CooChunkRange *ranges,
    int *col_size,
    int col_num,
    int chunk_num) {
    int max_chunk_size = splited_max_chunk_size(col_size, col_num, chunk_num);
    int i, current_chunk_size, current_chunk_begin, current_chunk_idx;
    for (i = 0,
        current_chunk_size = 0,
        current_chunk_begin = 0,
        current_chunk_idx = 0;
         i < col_num;
         ++i) {
        if (current_chunk_size + col_size[i] > max_chunk_size) {
            ranges[current_chunk_idx].row_begin = current_chunk_begin;
            ranges[current_chunk_idx].row_num = i - current_chunk_begin;
            ranges[current_chunk_idx].size = current_chunk_size;

            current_chunk_size = 0;
            current_chunk_begin = i;
            current_chunk_idx += 1;
        }

        current_chunk_size += col_size[i];
    }

    if (current_chunk_idx < chunk_num) {
        ranges[current_chunk_idx].row_begin = current_chunk_begin;
        ranges[current_chunk_idx].row_num = i - current_chunk_begin;
        ranges[current_chunk_idx].size = current_chunk_size;

        current_chunk_size = 0;
        current_chunk_begin = i;
        current_chunk_idx += 1;
    }

    while (current_chunk_idx < chunk_num) {
        ranges[current_chunk_idx].row_begin = col_num;
        ranges[current_chunk_idx].row_num = 0;
        ranges[current_chunk_idx].size = 0;
        current_chunk_idx += 1;
    }
}

struct Element {
    int row, column;
    double value;

    bool operator<(const Element &another) const {
        if (row != another.row) {
            return row < another.row;
        }

        return column < another.column;
    }
#ifdef CHECK_ELEMENT_DISPATCHING
    int chunk_idx, block_idx;
#endif
};

#ifdef USE_SIMD
    #include <immintrin.h>

void println256(const char *info, __m256i vec) {
    int *arr = (int *)&vec;
    fprintf(stderr, "%s ", info);
    for (int i = 0; i < 8; ++i) {
        fprintf(stderr, "%d ", arr[i]);
    }
    fprintf(stderr, "\n");
}

void simd_binary_search(
    int *res,
    int *col_begin,
    int *col_end,
    int size,
    int *keys) {
    __m256i one = _mm256_set1_epi32(1);
    __m256i all_ones = _mm256_set1_epi32(-1);

    __m256i key = _mm256_loadu_si256((__m256i *)keys);
    __m256i left = _mm256_set1_epi32(0);
    __m256i right = _mm256_set1_epi32(size);
    __m256i result = _mm256_set1_epi32(-1);
    __m256i mask3 = _mm256_set1_epi32(0);

    while (_mm256_movemask_epi8(mask3) != -1) {
        __m256i mid = _mm256_add_epi32(
            left,
            _mm256_srai_epi32(_mm256_sub_epi32(right, left), 1));
        __m256i begin = _mm256_i32gather_epi32(col_begin, mid, sizeof(int));
        __m256i end = _mm256_i32gather_epi32(col_end, mid, sizeof(int));
        __m256i mask1 =
            _mm256_xor_si256(_mm256_cmpgt_epi32(end, key), all_ones);
        __m256i mask2 = _mm256_cmpgt_epi32(begin, key);
        mask3 = _mm256_xor_si256(_mm256_or_si256(mask1, mask2), all_ones);
        left = _mm256_blendv_epi8(left, _mm256_add_epi32(mid, one), mask1);
        right = _mm256_blendv_epi8(right, mid, mask2);
        result = _mm256_blendv_epi8(result, mid, mask3);
    }

    _mm256_storeu_si256((__m256i *)res, result);
}

#else

void naive_binary_search(
    int *res,
    int *col_begin,
    int *col_end,
    int size,
    int *keys) {
    for (int i = 0; i < 8; ++i) {
        int l = 0, r = size;
        for (;;) {
            int m = l + (r - l) / 2;
            auto begin = col_begin[m];
            auto end = col_end[m];
            if (end <= keys[i]) {
                l = m + 1;
            } else if (keys[i] < begin) {
                r = m;
            } else {
                res[i] = m;
                break;
            }
        }
    }
}
#endif

struct {
    int num;
    int col_begin[64];
    int col_end[64];

    int data[99999];
} BSCache;

void init_bs_cache(void) {
    memset(BSCache.data, -1, sizeof(BSCache.data));
}

void cached_binary_search(
    int *res,
    int *col_begin,
    int *col_end,
    int size,
    int *keys) {
    for (int i = 0; i < 8; ++i) {
        if (BSCache.data[keys[i]] >= 0) {
            res[i] = BSCache.data[keys[i]];
        } else {
            int l = 0, r = size;
            for (;;) {
                int m = l + (r - l) / 2;
                auto begin = col_begin[m];
                auto end = col_end[m];
                if (end <= keys[i]) {
                    l = m + 1;
                } else if (keys[i] < begin) {
                    r = m;
                } else {
                    res[i] = m;
                    break;
                }
            }
            BSCache.data[keys[i]] = res[i];
        }
    }
}

int cached_binary_search_single(
    int *col_begin,
    int *col_end,
    int size,
    int key) {
    if (BSCache.data[key] >= 0) {
        return BSCache.data[key];
    } else {
        int l = 0, r = size;
        for (;;) {
            int m = l + (r - l) / 2;
            auto begin = col_begin[m];
            auto end = col_end[m];
            if (end <= key) {
                l = m + 1;
            } else if (key < begin) {
                r = m;
            } else {
                BSCache.data[key] = m;
                return m;
                break;
            }
        }
    }
}

void give_away_elements(
    const SplitedCooMatrix &result,
    int data_size,
    int *elements_row,
    int *elements_column,
    double *elements_value) {
    init_bs_cache();

    int *block_size =
        (int *)malloc(sizeof(int) * result.chunk_num * result.chunk_num);
    memset(block_size, 0, sizeof(int) * result.chunk_num * result.chunk_num);

    int *col_begin = (int *)malloc(sizeof(int) * result.chunk_num);
    int *col_end = (int *)malloc(sizeof(int) * result.chunk_num);
    for (int i = 0; i < result.chunk_num; ++i) {
        col_begin[i] = result.chunk_ranges[i].row_begin;
        col_end[i] = col_begin[i] + result.chunk_ranges[i].row_num;
    }

    int i;
    for (i = 0; i < data_size / 8 * 8; i += 8) {
        int block_idx[8];
        int chunk_idx[8];
#ifdef USE_SIMD
        simd_binary_search(
            chunk_idx,
            col_begin,
            col_end,
            result.chunk_num,
            &elements_row[i]);
        simd_binary_search(
            block_idx,
            col_begin,
            col_end,
            result.chunk_num,
            &elements_column[i]);
#else
        cached_binary_search(
            chunk_idx,
            col_begin,
            col_end,
            result.chunk_num,
            &elements_row[i]);
        cached_binary_search(
            block_idx,
            col_begin,
            col_end,
            result.chunk_num,
            &elements_column[i]);
#endif

        for (int j = 0; j < 8; ++j) {
            block_size[chunk_idx[j] * result.chunk_num + block_idx[j]] += 1;
        }
    }

    for (; i < data_size; ++i) {
        int chunk_idx = cached_binary_search_single(
            col_begin,
            col_end,
            result.chunk_num,
            elements_row[i]);
        int block_idx = cached_binary_search_single(
            col_begin,
            col_end,
            result.chunk_num,
            elements_column[i]);

        block_size[chunk_idx * result.chunk_num + block_idx] += 1;
    }

    int *chunk_block_off =
        (int *)malloc(sizeof(int) * result.chunk_num * result.chunk_num);
    memset(
        chunk_block_off,
        0,
        sizeof(int) * result.chunk_num * result.chunk_num);

    for (int chunk_idx = 0; chunk_idx < result.chunk_num; ++chunk_idx) {
        auto &chunk = *result.chunks[chunk_idx].chunk;
        int chunk_data_size = 0;
        for (int block_idx = 0; block_idx < result.chunk_num; ++block_idx) {
            auto &block = chunk.blocks[block_idx];
            int idx = chunk_idx * result.chunk_num + block_idx;
            block.block_off = chunk_data_size;
            chunk_block_off[idx] = chunk_data_size;
            chunk_data_size += block_size[idx];
        }
    }

    for (i = 0; i < data_size; ++i) {
        int block_idx = BSCache.data[elements_column[i]];
        int chunk_idx = BSCache.data[elements_row[i]];
        auto &chunk = *result.chunks[chunk_idx].chunk;
        auto &idx = chunk_block_off[chunk_idx * result.chunk_num + block_idx];

        chunk.row_idx[idx] = elements_row[i] - col_begin[chunk_idx];
        chunk.col_idx[idx] = elements_column[i] - col_begin[block_idx];
        chunk.data[idx] = elements_value[i];
        idx += 1;
    }

    free(chunk_block_off);
    free(block_size);
    free(col_begin);
    free(col_end);
}

extern void free_splited_coo(SplitedCooMatrix &mtx) {
    free(mtx.chunk_ranges);
    for (int i = 0; i < mtx.chunk_num; ++i) {
        auto &chunk = *mtx.chunks[i].chunk;
        free(chunk.col_idx);
        free(chunk.row_idx);
        free(chunk.data);
        free(mtx.chunks[i].chunk);
    }
    free(mtx.chunks);
}

extern SplitedCooMatrix
ldu_to_splited_coo(const LduMatrix &ldu_matrix, int chunk_num) {
    SplitedCooMatrix result;

    int row_num = ldu_matrix.cells;

    int *elements_row =
        (int *)malloc(sizeof(int) * (ldu_matrix.cells + ldu_matrix.faces * 2));
    int *elements_column =
        (int *)malloc(sizeof(int) * (ldu_matrix.cells + ldu_matrix.faces * 2));
    double *elements_value = (double *)malloc(
        sizeof(double) * (ldu_matrix.cells + ldu_matrix.faces * 2));
    int data_size = 0;

    int *row_size = (int *)malloc(sizeof(int) * row_num);

    for (int i = 0; i < ldu_matrix.cells; i++) {
        row_size[i] = 1;
    }
    for (int i = 0; i < ldu_matrix.faces; i++) {
        row_size[ldu_matrix.lPtr[i]] += 1;
        row_size[ldu_matrix.uPtr[i]] += 1;
    }

    memcpy(elements_row, ldu_matrix.uPtr, sizeof(int) * ldu_matrix.faces);
    memcpy(
        elements_row + ldu_matrix.faces,
        ldu_matrix.lPtr,
        sizeof(int) * ldu_matrix.faces);

    memcpy(elements_column, ldu_matrix.lPtr, sizeof(int) * ldu_matrix.faces);
    memcpy(
        elements_column + ldu_matrix.faces,
        ldu_matrix.uPtr,
        sizeof(int) * ldu_matrix.faces);

    memcpy(elements_value, ldu_matrix.lower, sizeof(double) * ldu_matrix.faces);
    memcpy(
        elements_value + ldu_matrix.faces,
        ldu_matrix.upper,
        sizeof(double) * ldu_matrix.faces);

    data_size = 2 * ldu_matrix.faces;

    memcpy(
        elements_value + data_size,
        ldu_matrix.diag,
        sizeof(double) * ldu_matrix.cells);

    for (int i = 0; i < ldu_matrix.cells; i++) {
        elements_row[data_size + i] = i;
        elements_column[data_size + i] = i;
    }

    data_size += ldu_matrix.cells;

    result.chunk_num = chunk_num;

    result.chunk_ranges =
        (CooChunkRange *)malloc(sizeof(CooChunkRange) * chunk_num);
    get_chunk_ranges(result.chunk_ranges, row_size, row_num, chunk_num);

    result.chunks = (SizedCooChunk *)malloc(sizeof(SizedCooChunk) * chunk_num);
    for (int i = 0; i < chunk_num; ++i) {
        result.chunks[i].mem_size =
            sizeof(CooChunk) + sizeof(CooBlock) * (chunk_num + 1);
        result.chunks[i].chunk = (CooChunk *)malloc(result.chunks[i].mem_size);

        auto &chunk = *result.chunks[i].chunk;

        int chunk_data_size = 0;
        for (int j = 0; j < result.chunk_ranges[i].row_num; ++j) {
            chunk_data_size += row_size[j + result.chunk_ranges[i].row_begin];
        }

        chunk.row_begin = result.chunk_ranges[i].row_begin;
        chunk.row_end = chunk.row_begin + result.chunk_ranges[i].row_num;
        chunk.row_idx = (uint16_t *)malloc(sizeof(uint16_t) * chunk_data_size);
        chunk.col_idx = (uint16_t *)malloc(sizeof(uint16_t) * chunk_data_size);
        chunk.data = (double *)malloc(sizeof(double) * chunk_data_size);
        chunk.block_num = chunk_num;

        for (int j = 0; j < chunk.block_num; ++j) {
            auto &block = chunk.blocks[j];
            block.row_num = chunk.row_end - chunk.row_begin;
            block.col_begin = result.chunk_ranges[j].row_begin;
            block.block_off = 0;
        }
        chunk.blocks[chunk.block_num].row_num = 0;
        chunk.blocks[chunk.block_num].col_begin = row_num;
        chunk.blocks[chunk.block_num].block_off = chunk_data_size;
    }

    give_away_elements(
        result,
        data_size,
        elements_row,
        elements_column,
        elements_value);

    free(row_size);
    free(elements_row);
    free(elements_column);
    free(elements_value);
    return result;
}
