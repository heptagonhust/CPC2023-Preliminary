#include "spmv_slave.h"
#include "crts.h"
#include "swperf.h"
#include "spmv_def.h"

#define SLAVE_CORE_NUM 64

__thread_local crts_rply_t rma_get_rply = 0;
__thread_local unsigned int rma_get_cnt = 0;
__thread_local crts_rply_t rma_remote_get_rply = 0;
__thread_local unsigned int rma_remote_get_cnt = 0;

inline void slave_double_buffering_new(DoubleBuffering *buff, int size) {
    buff->size = size;
    buff->buff[0] = (double *)CRTS_pldm_malloc(size * sizeof(double));
    buff->buff[1] = (double *)CRTS_pldm_malloc(size * sizeof(double));
    buff->in_use = 1;
}

inline void* slave_double_buffering_get(DoubleBuffering *buff) {
    int should_use = !buff->in_use;
    buff->in_use = should_use;
    return buff->buff[should_use];
}

inline void slave_double_buffering_free(DoubleBuffering *buff) {
    CRTS_pldm_free(buff->buff[0], buff->size * sizeof(double));
    CRTS_pldm_free(buff->buff[1], buff->size * sizeof(double));
}

inline void slave_coo_spmv(CooChunk *chunk, double *vec, double **vec_list, double *result, DoubleBuffering *buff, uint16_t non_0_block_num, uint16_t *non_0_block_idx) {
    int id = CRTS_tid;
    int rows = chunk->row_end - chunk->row_begin;
    memset(result, 0, rows * sizeof(double));

    int block_data_size, next_vec_size;
    double *vec_buf, *rma_buf;
    if (non_0_block_idx[0] == id) {
        block_data_size = chunk->blocks[id+1].block_off - chunk->blocks[id].block_off;
        if (non_0_block_num > 1) {
            next_vec_size = chunk->blocks[non_0_block_idx[1]+1].col_begin - chunk->blocks[non_0_block_idx[1]].col_begin;
            rma_buf = slave_double_buffering_get(buff);
            CRTS_rma_iget(rma_buf, &rma_get_rply, next_vec_size * sizeof(double), non_0_block_idx[1], vec_list[non_0_block_idx[1]], &rma_remote_get_rply);
            ++rma_get_cnt;
        }
        int block_data_offset = chunk->blocks[id].block_off;
        for (int j = 0; j < block_data_size; ++j) {
            result[chunk->row_idx[block_data_offset + j]] += chunk->data[block_data_offset + j] * vec[chunk->col_idx[block_data_offset + j]];
        }
    }

    for (int i = 1; i < non_0_block_num; ++i) {
        CRTS_rma_wait_value(&rma_get_rply, rma_get_cnt);
        vec_buf = rma_buf;
        block_data_size = chunk->blocks[non_0_block_idx[i] + 1].block_off - chunk->blocks[non_0_block_idx[i]].block_off;
        if (i < non_0_block_num - 1) {
            next_vec_size = chunk->blocks[non_0_block_idx[i+1] + 1].col_begin - chunk->blocks[non_0_block_idx[i+1]].col_begin; 
            rma_buf = slave_double_buffering_get(buff);
            CRTS_rma_iget(rma_buf, &rma_get_rply, next_vec_size * sizeof(double), non_0_block_idx[i+1], vec_list[non_0_block_idx[i+1]], &rma_remote_get_rply);
            ++rma_get_cnt;
        }

        int block_data_offset = chunk->blocks[non_0_block_idx[i]].block_off;
        for (int j = 0; j < block_data_size; ++j) {
            result[chunk->row_idx[block_data_offset + j]] += chunk->data[block_data_offset + j] * vec_buf[chunk->col_idx[block_data_offset + j]];
        }
    }
    CRTS_ssync_array();
}