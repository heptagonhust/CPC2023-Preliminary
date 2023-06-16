#include "coo_matrix.h"

#include <algorithm>
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
            ranges[current_chunk_idx].col_begin = current_chunk_begin;
            ranges[current_chunk_idx].col_num = i - current_chunk_begin;
            ranges[current_chunk_idx].size = current_chunk_size;

            current_chunk_size = 0;
            current_chunk_begin = i;
            current_chunk_idx += 1;
        }

        current_chunk_size += col_size[i];
    }

    if (current_chunk_idx < chunk_num) {
        ranges[current_chunk_idx].col_begin = current_chunk_begin;
        ranges[current_chunk_idx].col_num = i - current_chunk_begin;
        ranges[current_chunk_idx].size = current_chunk_size;

        current_chunk_size = 0;
        current_chunk_begin = i;
        current_chunk_idx += 1;
    }

    while (current_chunk_idx < chunk_num) {
        ranges[current_chunk_idx].col_begin = col_num;
        ranges[current_chunk_idx].col_num = 0;
        ranges[current_chunk_idx].size = 0;
        current_chunk_idx += 1;
    }
}

extern SplitedCooMatrix
ldu_to_splited_coo(const LduMatrix &ldu_matrix, int chunk_num) {
    SplitedCooMatrix result;

    int col_num = ldu_matrix.cells;

    typedef struct Element {
        int row, column;
        int chunk_idx, block_idx;
        double value;

        bool operator<(const Element &another) {
            if (chunk_idx != another.chunk_idx) {
                return chunk_idx < another.chunk_idx;
            }

            if (block_idx != another.block_idx) {
                return block_idx < another.block_idx;
            }

            if (row != another.row) {
                return row < another.row;
            }

            return column < another.column;
        }
    } Element;

    std::vector<Element> elements;
    elements.reserve(col_num / 20 * col_num);

    int *col_size = (int *)malloc(sizeof(int) * col_num);
    memset(col_size, 0, sizeof(int) * col_num);

    for (int i = 0; i < ldu_matrix.faces; i++) {
        int row, column;
        double value;

        value = ldu_matrix.lower[i];
        if (value != 0.0) {
            row = ldu_matrix.uPtr[i];
            column = ldu_matrix.lPtr[i];
            elements.push_back((Element) {row, column, -1, -1, value});
            col_size[column] += 1;
        }

        value = ldu_matrix.upper[i];
        if (value != 0.0) {
            row = ldu_matrix.lPtr[i];
            column = ldu_matrix.uPtr[i];
            elements.push_back((Element) {row, column, -1, -1, value});
            col_size[column] += 1;
        }
    }

    for (int i = 0; i < ldu_matrix.cells; i++) {
        auto value = ldu_matrix.diag[i];
        if (value != 0.0) {
            auto row = i;
            auto column = i;
            elements.push_back((Element) {row, column, -1, -1, value});
            col_size[column] += 1;
        }
    }

    result.chunk_num = chunk_num;

    result.chunk_ranges =
        (CooChunkRange *)malloc(sizeof(CooChunkRange) * chunk_num);
    get_chunk_ranges(result.chunk_ranges, col_size, col_num, chunk_num);

    result.chunks = (SizedCooChunk *)malloc(sizeof(SizedCooChunk) * chunk_num);
    for (int i = 0; i < chunk_num; ++i) {
        result.chunks[i].mem_size =
            sizeof(CooChunk) + sizeof(CooBlock) * (chunk_num + 1);
        result.chunks[i].chunk = (CooChunk *)malloc(result.chunks[i].mem_size);

        auto &chunk = *result.chunks[i].chunk;

        int data_size = 0;
        for (int j = 0; j < result.chunk_ranges[i].col_num; ++j) {
            data_size += col_size[j + result.chunk_ranges[i].col_begin];
        }

        chunk.col_begin = result.chunk_ranges[i].col_begin;
        chunk.col_end = chunk.col_begin + result.chunk_ranges[i].col_num;
        chunk.row_idx = (int *)malloc(sizeof(int) * data_size);
        chunk.col_idx = (int *)malloc(sizeof(int) * data_size);
        chunk.data = (double *)malloc(sizeof(double) * data_size);
        chunk.block_num = chunk_num;

        for (int j = 0; j < chunk.block_num; ++j) {
            auto &block = chunk.blocks[j];
            block.col_num = chunk.col_end - chunk.col_begin;
            block.row_begin = result.chunk_ranges[j].col_begin;
            block.block_off = 0;
        }
        chunk.blocks[chunk.block_num].col_num = 0;
        chunk.blocks[chunk.block_num].row_begin = col_num;
        chunk.blocks[chunk.block_num].block_off = data_size;
    }

    for (auto element = elements.begin(); element != elements.end();
         ++element) {
        int chunk_idx = -1;
        int block_idx = -1;
        for (int i = 0; i < chunk_num; ++i) {
            auto begin = result.chunk_ranges[i].col_begin;
            auto end = begin + result.chunk_ranges[i].col_num;
            if (begin <= element->row && element->row < end) {
                block_idx = i;
            }
            if (begin <= element->column && element->column < end) {
                chunk_idx = i;
            }
            if (block_idx >= 0 && chunk_idx >= 0) {
                break;
            }
        }
        element->chunk_idx = chunk_idx;
        element->block_idx = block_idx;
    }

    std::sort(elements.begin(), elements.end());

    int chunk_data_size = 0, last_chunk_idx = -1;
    for (auto element = elements.begin(); element != elements.end();
         ++element) {
        if (element->chunk_idx != last_chunk_idx) {
            chunk_data_size = 0;
        }

        auto &chunk = *result.chunks[element->chunk_idx].chunk;
        auto &block = chunk.blocks[element->block_idx];
        chunk.row_idx[chunk_data_size] = element->row - block.row_begin;
        chunk.col_idx[chunk_data_size] = element->column - chunk.col_begin;
        chunk.data[chunk_data_size] = element->value;
        chunk_data_size += 1;
        block.block_off += 1;

        last_chunk_idx = element->chunk_idx;
    }

    for (int i = 0; i < chunk_num; ++i) {
        auto &chunk = *result.chunks[i].chunk;
        for (int j = 1; j < chunk.block_num; ++j) {
            chunk.blocks[j].block_off += chunk.blocks[j - 1].block_off;
        }
        for (int j = chunk.block_num - 1; j >= 1; --j) {
            chunk.blocks[j].block_off -=
                (chunk.blocks[j].block_off - chunk.blocks[j - 1].block_off);
        }
        chunk.blocks[0].block_off = 0;
    }

    return result;
}