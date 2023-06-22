#include "spmv_slave.h"
#include "crts.h"
#include "swperf.h"
#include "spmv_def.h"


inline void slave_coo_spmv(CooChunk *chunk, double *vec, double *result) {
    int rows = chunk->row_end - chunk->row_begin;
    memset(result, 0, rows * sizeof(double));

    intv16 off_a, off_b;
    for (int i = 0; i < chunk->block_num; ++i) {
        if (i % 16 == 0) {
            simd_loadu(off_b, chunk->block_off + i);
            simd_loadu(off_a, chunk->block_off + i + 1);
            off_a -= off_b;
        }
        int block_data_size = ((int *)(&off_a))[i % 16];
        if (block_data_size != 0) {
            int block_data_offset = chunk->block_off[i];
            int block_col_begin = chunk->blocks[i].col_begin;
            for (int j = 0; j < block_data_size; ++j) {
                result[chunk->row_idx[block_data_offset + j]] += chunk->data[block_data_offset + j] * vec[block_col_begin + chunk->col_idx[block_data_offset + j]];
            }
        }
    }
}