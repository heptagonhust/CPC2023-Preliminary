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
static int splited_max_chunk_size(
    int *col_size,
    int total_size,
    int col_num,
    int chunk_num) {
    int l = 0, r = total_size;
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

static bool real_get_chunk_ranges(
    CooChunkRange *ranges,
    int *col_size,
    int col_num,
    int chunk_num,
    int max_chunk_size) {
    int i, current_chunk_size, current_chunk_begin, current_chunk_idx;
    for (i = 0,
        current_chunk_size = 0,
        current_chunk_begin = 0,
        current_chunk_idx = 0;
         i < col_num;
         ++i) {
        if (current_chunk_idx >= chunk_num) {
            return false;
        }
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

    return true;
}

struct {
    int num;
    int col_begin[64];
    int col_end[64];

    int data[99999];
} BSCache;

void init_bs_cache(void) {
    memset(BSCache.data, -1, sizeof(BSCache.data));
}

static void get_chunk_ranges_and_init_bs_cache(
    CooChunkRange *ranges,
    int *col_size,
    int total_size,
    int col_num,
    int chunk_num) {
    static int last_max_chunk_size;

    if (last_max_chunk_size != 0) {
        if (real_get_chunk_ranges(
                ranges,
                col_size,
                col_num,
                chunk_num,
                last_max_chunk_size)) {
            int max_chunk = 0;
            int min_chunk = 999999999;
            for (int i = 0; i < chunk_num; ++i) {
                min_chunk = std::min(min_chunk, ranges[i].size);
                max_chunk = std::max(max_chunk, ranges[i].size);
            }
            if (8 * (max_chunk - min_chunk) < max_chunk) {
                return;
            }
        }
    }

    last_max_chunk_size =
        splited_max_chunk_size(col_size, total_size, col_num, chunk_num);
    real_get_chunk_ranges(
        ranges,
        col_size,
        col_num,
        chunk_num,
        last_max_chunk_size);
    init_bs_cache();
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

    int *mem = (int *)malloc(std::max(
        sizeof(int) * ldu_matrix.cells,
        sizeof(int) * (chunk_num * 2 + chunk_num * chunk_num * 2)));
    int *col_begin = mem;
    int *col_end = mem + chunk_num;
    int *chunk_block_off = mem + chunk_num * 2;
    int *block_size = mem + chunk_num * 2 + chunk_num * chunk_num;
    int *row_size = mem;

    result.chunk_num = chunk_num;
    result.chunk_ranges =
        (CooChunkRange *)malloc(sizeof(CooChunkRange) * chunk_num);

    int total_size = 0;
    for (int i = 0; i < ldu_matrix.cells; i++) {
        row_size[i] = 1;
    }
    for (int i = 0; i < ldu_matrix.faces; i++) {
        row_size[ldu_matrix.lPtr[i]] += 1;
        row_size[ldu_matrix.uPtr[i]] += 1;
    }

    for (int i = 0; i < row_num; ++i) {
        total_size += row_size[i];
    }

    get_chunk_ranges_and_init_bs_cache(
        result.chunk_ranges,
        row_size,
        total_size,
        row_num,
        chunk_num);

    result.chunks = (SizedCooChunk *)malloc(sizeof(SizedCooChunk) * chunk_num);
    for (int i = 0; i < chunk_num; ++i) {
        result.chunks[i].mem_size =
            sizeof(CooChunk) + sizeof(CooBlock) * (chunk_num + 1);
        result.chunks[i].chunk = (CooChunk *)malloc(result.chunks[i].mem_size);

        auto &chunk = *result.chunks[i].chunk;

        int chunk_data_size = result.chunk_ranges[i].size;

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

    for (int i = 0; i < result.chunk_num; ++i) {
        col_begin[i] = result.chunk_ranges[i].row_begin;
        col_end[i] = col_begin[i] + result.chunk_ranges[i].row_num;
    }

    memset(block_size, 0, sizeof(int) * chunk_num * chunk_num);
    for (int i = 0; i < ldu_matrix.faces; ++i) {
        int chunk_idx = cached_binary_search_single(
            col_begin,
            col_end,
            result.chunk_num,
            ldu_matrix.uPtr[i]);
        int block_idx = cached_binary_search_single(
            col_begin,
            col_end,
            result.chunk_num,
            ldu_matrix.lPtr[i]);
        block_size[chunk_idx * result.chunk_num + block_idx] += 1;
        block_size[block_idx * result.chunk_num + chunk_idx] += 1;
    }

    for (int i = 0; i < ldu_matrix.cells; ++i) {
        int idx = cached_binary_search_single(
            col_begin,
            col_end,
            result.chunk_num,
            i);
        block_size[idx * result.chunk_num + idx] += 1;
    }

    memset(chunk_block_off, 0, sizeof(int) * chunk_num * chunk_num);
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

    for (int i = 0; i < ldu_matrix.faces; ++i) {
        int idx1 = BSCache.data[ldu_matrix.uPtr[i]];
        int idx2 = BSCache.data[ldu_matrix.lPtr[i]];
        int idx_21 = chunk_block_off[idx2 * result.chunk_num + idx1];
        chunk_block_off[idx2 * result.chunk_num + idx1] += 1;
        int idx_12 = chunk_block_off[idx1 * result.chunk_num + idx2];
        chunk_block_off[idx1 * result.chunk_num + idx2] += 1;

        result.chunks[idx2].chunk->row_idx[idx_21] =
            ldu_matrix.lPtr[i] - col_begin[idx2];
        result.chunks[idx2].chunk->col_idx[idx_21] =
            ldu_matrix.uPtr[i] - col_begin[idx1];
        result.chunks[idx2].chunk->data[idx_21] = ldu_matrix.upper[i];

        result.chunks[idx1].chunk->row_idx[idx_12] =
            ldu_matrix.uPtr[i] - col_begin[idx1];
        result.chunks[idx1].chunk->col_idx[idx_12] =
            ldu_matrix.lPtr[i] - col_begin[idx2];
        result.chunks[idx1].chunk->data[idx_12] = ldu_matrix.lower[i];
    }

    for (int i = 0; i < ldu_matrix.cells; ++i) {
        int id = BSCache.data[i];
        auto &chunk = *result.chunks[id].chunk;
        auto &idx = chunk_block_off[id * result.chunk_num + id];
        chunk.row_idx[idx] = i - col_begin[id];
        chunk.col_idx[idx] = i - col_begin[id];
        chunk.data[idx] = ldu_matrix.diag[i];
        idx += 1;
    }

    free(mem);

    return result;
}
