#include "spmv_master.h"

#define SLAVE_CORE_NUM 64
void csc_spmv(SpmvPara *para, double *result) {
    CRTS_athread_spawn(slave_csc_spmv, &para);
    int chunk_num = para->chunk_num;
    double **para_result = para->result;
    int *dma_over = para->dma_over;
    CscChunk *chunk0 = para->chunks[0];
    int block_num_per_chunk = chunk0->block_num;

    for (int i = 0; i < block_num_per_chunk; ++i) {
        for (int j = 0; j < chunk_num; ++j) {
            while (dma_over[j] < i) {
                // wait for dma
            }
        }
        for (int j = 0; j < chunk_num; ++j) {
            CscBlock *block = chunk0->blocks[i];
            int row_begin = block->row_begin;
            int row_num = block->row_end - row_begin;
            double *chunk_result = para_result[j] + row_begin;
            result += row_begin;
            for (int k = 0; k < row_num; ++k) {
                result[k] += chunk_result[k];
            }
        }
    }
}