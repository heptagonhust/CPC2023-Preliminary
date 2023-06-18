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
    int pool_offset = 0;
    int reduced_cnt = 0;
    bool reduce_status[chunk_num * block_num_per_chunk];
    memset(&reduce_status, 0, chunk_num * block_num_per_chunk * sizeof(bool));

    while (reduced_cnt < chunk_num * block_num_per_chunk) {
        for (int block_idx = 0; block_idx < block_num_per_chunk; ++i) {
            for (int chunk_idx = 0; chunk_idx < chunk_num; ++j) {
                int flag_idx = chunk_num * block_idx + chunk_idx;
                if (dma_over[chunk_idx] >= block_idx && !reduce_status[flag_idx]) {
                    reduce_status[flag_idx] = 1;
                    ++reduced_cnt;
                }
            }
        }
    }

    for (int i = 0; i < block_num_per_chunk; ++i) {
        // 每次规约前要等待所有从核 DMA 完成
        CscBlock *block = chunk0->blocks + i;
        int row_begin = block->row_begin;
        int row_num = block->row_end - row_begin;

        reduced_cnt = 0;
        
        while (reduced_cnt < chunk_num) {
            for (int j = 0; j < chunk_num; ++j) {
                if (dma_over[j] >= i && !reduce_status[j]) {
                    // wait for dma
                    for (int k = 0; k < row_num; ++k, ++result) {
                        result[row_begin + k] += result_pool[pool_offset + k * chunk_num + j];
                    }
                    reduce_status[j] = 1;
                    ++reduced_cnt;
                }
            }
        }
        pool_offset += row_num * chunk_num;
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