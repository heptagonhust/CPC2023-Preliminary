#include "spmv_master.h"

#define SLAVE_CORE_NUM 64
void csc_spmv(SpmvPara *para, double *result) {
    CRTS_athread_spawn(slave_csc_spmv, para);
    memset(result, 0, sizeof(double) * para->sp_col);
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
        CscBlock *block = chunk0->blocks + i;
        int row_begin = block->row_begin;
        int row_num = block->row_end - row_begin;
        for (int j = 0; j < row_num; ++j, ++result) {
            for (int k = 0; k < chunk_num; ++k, ++result_pool) {
                *result += *result_pool;
            }
        }
    }
}

void spmv_para_from_splited_coo_matrix(
    const SplitedCscMatrix *mat,
    SpmvPara *para,
    double *vec,
    int row_num,
    int col_num) {
    int chunk_num = mat->chunk_num;
    para->chunk_num = chunk_num;
    para->sp_row = row_num;
    para->sp_col = col_num;
    para->result = (double *)malloc(sizeof(double) * chunk_num * row_num);
    para->dma_over = (int *)malloc(sizeof(int) * chunk_num);
    memset(para->dma_over, 0, sizeof(int) * chunk_num);
    para->vec = vec;
    memcpy(para->chunks, mat->chunks, SLAVE_CORE_NUM * sizeof(SizedCooChunk));
    // all chunks have the same block size
    CooChunk *chunk0 = mat->chunks[0];
    int max_block_row_num = 0;
    for (int i = 0; i < chunk0->block_num; ++i) {
        int num = chunk0->blocks[i+1].row_begin - chunk0->blocks[i].row_begin;
        if (num > max_block_row_num) {
            max_block_row_num = num;
        }
    }
    para->max_block_row_num = max_block_row_num;
}

void spmv_para_free(SpmvPara *para) {
    free(para->dma_over);
    free(para->result);
}