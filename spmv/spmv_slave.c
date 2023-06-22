#include "spmv_slave.h"
#include "crts.h"
#include "swperf.h"
#include "spmv_def.h"

extern __thread_local unsigned long icc;
extern __thread_local unsigned long share_access_cycles;
extern __thread_local unsigned long pldm_access_cycles;


inline void slave_coo_spmv(CooChunk *chunk, double *vec, double *result) {
    int rows = chunk->row_end - chunk->row_begin;
    memset(result, 0, rows * sizeof(double));

    for (int i = 0; i < chunk->block_num; ++i) {
        int block_data_size = chunk->blocks[i+1].block_off - chunk->blocks[i].block_off;
        if (block_data_size != 0) {
            int block_data_offset = chunk->blocks[i].block_off;
            int block_col_begin = chunk->blocks[i].col_begin;
            for (int j = 0; j < block_data_size; ++j) {
                penv_slave0_cycle_count(&icc);
                int start = icc;
                double val_s = vec[block_col_begin + chunk->col_idx[block_data_offset + j]];
                penv_slave0_cycle_count(&icc);
                share_access_cycles += icc - start;

                penv_slave0_cycle_count(&icc);
                start = icc;
                double val = chunk->data[block_data_offset + j];
                penv_slave0_cycle_count(&icc);
                pldm_access_cycles += icc - start;
                result[chunk->row_idx[block_data_offset + j]] += val * val_s;
            }
        }
    }
}