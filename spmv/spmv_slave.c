#include "spmv_slave.h"
#include "crts.h"
#include "swperf.h"
#include "spmv_def.h"


inline void slave_coo_spmv(CooChunk *chunk, double *vec, double *result) {
    int rows = chunk->row_end - chunk->row_begin;
    memset(result, 0, rows * sizeof(double));

    for (int i = 0; i < chunk->block_num; ++i) {
        int block_data_size = chunk->blocks[i+1].block_off - chunk->blocks[i].block_off;
        if (block_data_size != 0) {
            int block_data_offset = chunk->blocks[i].block_off;
            int block_col_begin = chunk->blocks[i].col_begin;
            for (int j = 0; j < block_data_size; ++j) {
                result[chunk->row_idx[block_data_offset + j]] += chunk->data[block_data_offset + j] * vec[block_col_begin + chunk->col_idx[block_data_offset + j]];
            }
        }
    }
}