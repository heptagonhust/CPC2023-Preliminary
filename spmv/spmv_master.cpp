#include "spmv_master.h"

#define SLAVE_CORE_NUM 64
void csc_spmv(SpmvPara *para, double *result) {
    CRTS_athread_spawn(slave_csc_spmv, &para);
    int chunk_num = para->chunk_num;
    double *result_pool = para->result;
    int *dma_over = para->dma_over;
    CscChunk *chunk0 = para->chunks[0];
    int block_num_per_chunk = chunk0->block_num;

    for (int i = 0; i < block_num_per_chunk; ++i) {
        // 每次规约前要等待所有从核 DMA 完成
        for (int j = 0; j < chunk_num; ++j) {
            while (dma_over[j] < i) {
                // wait for dma
            }
        }
        // 开始规约
        CscBlock *block = chunk0->blocks[i];
        int row_begin = block->row_begin;
        int row_num = block->row_end - row_begin;
        for (int j = 0; j < row_num; ++j) {
            for (int k = 0; k < chunk_num; ++k, ++result_pool) {
                *result += *result_pool;
            }
            result += 1;
        }
    }
}