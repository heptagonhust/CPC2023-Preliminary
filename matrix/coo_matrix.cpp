#include "coo_matrix.h"

#include <math.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "matrix_utils.h"
#include "pcg.h"
#include "spmv_def.h"

#define SHADOW_BLOCK_SIZE 256

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
        int chunk_cnt = 0;
        for (int i = 0, current_chunk_size = 0; i < col_num; ++i) {
            current_chunk_size += col_size[i];
            if (current_chunk_size > chunk_size) {
                if (col_size[i] > chunk_size) {
                    chunk_cnt = 0x3f3f3f3f;
                    break;
                }
                current_chunk_size = col_size[i];
                chunk_cnt += 1;
            }
        }
        chunk_cnt += 1;

        if (chunk_cnt <= chunk_num) {
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
    int *col_shadow,
    int chunk_num,
    int max_chunk_size) {
    int i = 0, current_chunk_size = 0, current_chunk_begin = 0,
        current_chunk_idx = 0;

    while (i < col_num) {
        if (i + SHADOW_BLOCK_SIZE < col_num) {
            int new_chunk_size =
                current_chunk_size + col_shadow[i / SHADOW_BLOCK_SIZE];

            if (new_chunk_size <= max_chunk_size) {
                current_chunk_size = new_chunk_size;
            } else {
                for (int j = 0; j < SHADOW_BLOCK_SIZE; ++j) {
                    if (col_size[i + j] > max_chunk_size) {
                        return false;
                    }

                    if (current_chunk_size + col_size[i + j] > max_chunk_size) {
                        if (current_chunk_idx >= chunk_num) {
                            return false;
                        }
                        ranges[current_chunk_idx].row_begin =
                            current_chunk_begin;
                        ranges[current_chunk_idx].row_num =
                            i + j - current_chunk_begin;
                        ranges[current_chunk_idx].size = current_chunk_size;

                        current_chunk_size = 0;
                        current_chunk_begin = i + j;
                        current_chunk_idx += 1;
                    }
                    current_chunk_size += col_size[i + j];
                }
            }

            i += SHADOW_BLOCK_SIZE;
        } else {
            if (col_size[i] > max_chunk_size) {
                return false;
            }

            if (current_chunk_size + col_size[i] > max_chunk_size) {
                ranges[current_chunk_idx].row_begin = current_chunk_begin;
                ranges[current_chunk_idx].row_num = i - current_chunk_begin;
                ranges[current_chunk_idx].size = current_chunk_size;

                current_chunk_size = 0;
                current_chunk_begin = i;
                if (current_chunk_idx >= chunk_num) {
                    return false;
                }
                current_chunk_idx += 1;
            }

            current_chunk_size += col_size[i];
            i += 1;
        }
    }

    if (current_chunk_size != 0) {
        ranges[current_chunk_idx].row_begin = current_chunk_begin;
        ranges[current_chunk_idx].row_num = i - current_chunk_begin;
        ranges[current_chunk_idx].size = current_chunk_size;

        current_chunk_size = 0;
        current_chunk_begin = i;
        current_chunk_idx += 1;
    }

    if (current_chunk_idx > chunk_num) {
        return false;
    }

    while (current_chunk_idx < chunk_num) {
        ranges[current_chunk_idx].row_begin = col_num;
        ranges[current_chunk_idx].row_num = 0;
        ranges[current_chunk_idx].size = 0;
        current_chunk_idx += 1;
    }

    return true;
}

static void get_chunk_ranges(
    CooChunkRange *ranges,
    int *col_size,
    int total_size,
    int col_num,
    int *col_shadow,
    int chunk_num) {
    int max_chunk_size =
        splited_max_chunk_size(col_size, total_size, col_num, chunk_num);

    real_get_chunk_ranges(
        ranges,
        col_size,
        col_num,
        col_shadow,
        chunk_num,
        max_chunk_size);
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

    int *mem = (int *)malloc(sizeof(int) * (chunk_num * chunk_num * 2));
    int *chunk_block_off = mem;
    int *block_size = chunk_block_off + chunk_num * chunk_num;

    CooChunk **chunks = (CooChunk **)malloc(sizeof(CooChunk *) * chunk_num);
    CooBlock **blocks =
        (CooBlock **)malloc(sizeof(CooBlock *) * chunk_num * chunk_num);
    CooChunkRange **ranges =
        (CooChunkRange **)malloc(sizeof(CooChunkRange *) * chunk_num);
    int *uPtr = ldu_matrix.uPtr;
    int *lPtr = ldu_matrix.lPtr;
    double *diag = ldu_matrix.diag;

    static int bscache[99999];
    static int last_matrix_size;
    static CooChunkRange last_split[99];

    int hit = true;

    if (last_matrix_size != row_num) {
        hit = false;
    }

    if (hit) {
        int *chunk_size = (int *)malloc(sizeof(int) * chunk_num);
        memset(chunk_size, 0, sizeof(int) * chunk_num);
        for (int i = 0; i < ldu_matrix.faces; i += 512) {
            chunk_size[bscache[uPtr[i]]] += 1;
            chunk_size[bscache[lPtr[i]]] += 1;
        }
        int max_chunk = 0;
        int min_chunk = 999999999;
        for (int i = 0; i < chunk_num; ++i) {
            min_chunk = std::min(min_chunk, chunk_size[i]);
            max_chunk = std::max(max_chunk, chunk_size[i]);
        }
        if (10 * (max_chunk - min_chunk) < max_chunk) {
            hit = false;
        }
        free(chunk_size);
    }

    result.chunk_num = chunk_num;
    result.chunk_ranges =
        (CooChunkRange *)malloc(sizeof(CooChunkRange) * chunk_num);

    if (hit) {
        memcpy(
            result.chunk_ranges,
            last_split,
            sizeof(CooChunkRange) * chunk_num);
    } else {
        int *tmp_mem = (int *)malloc(
            sizeof(int) * (ldu_matrix.cells + row_num / SHADOW_BLOCK_SIZE + 1));
        int *row_size = tmp_mem;
        int *row_shadow = tmp_mem + ldu_matrix.cells;

        memset(row_shadow, 0, sizeof(int) * (row_num / SHADOW_BLOCK_SIZE + 1));
        memset(row_size, 0, sizeof(int) * row_num);
        int total_size = 0;
        for (int i = 0; i < ldu_matrix.cells; i++) {
            row_size[i] = 1;
            row_shadow[i / SHADOW_BLOCK_SIZE] += 1;
        }
        for (int i = 0; i < ldu_matrix.faces; i++) {
            row_size[lPtr[i]] += 1;
            row_size[uPtr[i]] += 1;
            row_shadow[lPtr[i] / SHADOW_BLOCK_SIZE] += 1;
            row_shadow[uPtr[i] / SHADOW_BLOCK_SIZE] += 1;
        }

        for (int i = 0; i < row_num / SHADOW_BLOCK_SIZE + 1; ++i) {
            total_size += row_shadow[i];
        }

        get_chunk_ranges(
            result.chunk_ranges,
            row_size,
            total_size,
            row_num,
            row_shadow,
            chunk_num);

        for (int i = 0; i < chunk_num; ++i) {
            ranges[i] = &result.chunk_ranges[i];
            for (int j = ranges[i]->row_begin;
                 j < ranges[i]->row_begin + ranges[i]->row_num;
                 ++j) {
                bscache[j] = i;
            }
        }

        last_matrix_size = ldu_matrix.cells;
        memcpy(
            last_split,
            result.chunk_ranges,
            sizeof(CooChunkRange) * chunk_num);

        free(tmp_mem);
    }

    result.chunks = (SizedCooChunk *)malloc(sizeof(SizedCooChunk) * chunk_num);
    for (int i = 0; i < chunk_num; ++i) {
        result.chunks[i].mem_size =
            sizeof(CooChunk) + sizeof(CooBlock) * (chunk_num + 1);
        result.chunks[i].chunk = (CooChunk *)malloc(result.chunks[i].mem_size);

        chunks[i] = result.chunks[i].chunk;
        ranges[i] = &result.chunk_ranges[i];

        int chunk_data_size = ranges[i]->size;
        CooChunk *chunk = chunks[i];

        chunk->row_begin = result.chunk_ranges[i].row_begin;
        chunk->row_end = chunk->row_begin + result.chunk_ranges[i].row_num;
        chunk->row_idx = (uint16_t *)malloc(sizeof(uint16_t) * chunk_data_size);
        chunk->col_idx = (uint16_t *)malloc(sizeof(uint16_t) * chunk_data_size);
        chunk->data = (double *)malloc(sizeof(double) * chunk_data_size);
        chunk->block_num = chunk_num;

        for (int j = 0; j < chunk_num; ++j) {
            blocks[i * chunk_num + j] = &chunk->blocks[j];
            CooBlock *block = blocks[i * chunk_num + j];
            block->row_num = chunk->row_end - chunk->row_begin;
            block->col_begin = result.chunk_ranges[j].row_begin;
            block->block_off = 0;
        }
        chunk->blocks[chunk_num].row_num = 0;
        chunk->blocks[chunk_num].col_begin = row_num;
        chunk->blocks[chunk_num].block_off = chunk_data_size;
    }

    memset(block_size, 0, sizeof(int) * chunk_num * chunk_num);

    for (int i = 0; i < chunk_num; ++i) {
        block_size[i * chunk_num + i] = ranges[i]->row_num;
    }

    for (int i = 0; i < ldu_matrix.faces; ++i) {
        int chunk_idx = bscache[uPtr[i]];
        int block_idx = bscache[lPtr[i]];
        block_size[chunk_idx * chunk_num + block_idx] += 1;
        block_size[block_idx * chunk_num + chunk_idx] += 1;
    }

    memset(chunk_block_off, 0, sizeof(int) * chunk_num * chunk_num);
    for (int chunk_idx = 0; chunk_idx < chunk_num; ++chunk_idx) {
        auto &chunk = *chunks[chunk_idx];
        int chunk_data_size = 0;
        for (int block_idx = 0; block_idx < result.chunk_num; ++block_idx) {
            auto &block = chunk.blocks[block_idx];
            int idx = chunk_idx * result.chunk_num + block_idx;
            block.block_off = chunk_data_size;
            chunk_block_off[idx] = chunk_data_size;
            chunk_data_size += block_size[idx];
        }
    }

    for (int i = 0; i < chunk_num; ++i) {
        auto &idx = chunk_block_off[i * chunk_num + i];
        auto &chunk = *result.chunks[i].chunk;
        memcpy(
            chunks[i]->data + idx,
            diag + ranges[i]->row_begin,
            sizeof(double) * ranges[i]->row_num);
        for (int j = 0; j < ranges[i]->row_num; ++j) {
            chunk.row_idx[idx] = j;
            chunk.col_idx[idx] = j;
            idx += 1;
        }
    }

    for (int i = 0; i < ldu_matrix.faces; ++i) {
        int idx1 = bscache[uPtr[i]];
        int idx1_offset = uPtr[i] - ranges[idx1]->row_begin;
        int idx2 = bscache[lPtr[i]];
        int idx2_offset = lPtr[i] - ranges[idx2]->row_begin;
        int idx_21 = chunk_block_off[idx2 * chunk_num + idx1];
        chunk_block_off[idx2 * chunk_num + idx1] += 1;
        int idx_12 = chunk_block_off[idx1 * chunk_num + idx2];
        chunk_block_off[idx1 * chunk_num + idx2] += 1;

        chunks[idx1]->row_idx[idx_12] = idx1_offset;
        chunks[idx1]->col_idx[idx_12] = idx2_offset;
        chunks[idx1]->data[idx_12] = ldu_matrix.lower[i];

        chunks[idx2]->row_idx[idx_21] = idx2_offset;
        chunks[idx2]->col_idx[idx_21] = idx1_offset;
        chunks[idx2]->data[idx_21] = ldu_matrix.upper[i];
    }

    free(mem);
    free(chunks);
    free(blocks);
    free(ranges);

    return result;
}
