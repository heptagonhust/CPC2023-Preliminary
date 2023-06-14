#include "spmv_master.h"
#include <sys/_types/_int64_t.h>
#include <sys/_types/_size_t.h>
#include <sys/_types/_u_int64_t.h>
#include <cstring>

#define SLAVE_CORE_NUM 64
void csc_spmv(SpmvPara *para, double *result) {
    CRTS_athread_spawn(slave_csc_spmv, &para);
    int chunk_num = para->chunk_num;
    double *result_pool = para->result;
    int *dma_over = para->dma_over;
    CscChunk *chunk0 = para->chunks[0];
    int block_num_per_chunk = chunk0->block_num;
    int reduced_cnt = 0;
    int pool_offset = 0;
    bool reduce_status[chunk_num];

    for (int i = 0; i < block_num_per_chunk; ++i) {
        // 每次规约前要等待所有从核 DMA 完成
        CscBlock *block = chunk0->blocks + i;
        int row_begin = block->row_begin;
        int row_num = block->row_end - row_begin;

        reduced_cnt = 0;
        memset(&reduce_status, 0, chunk_num * sizeof(bool));
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

void splited_csc_matrix_to_spmv_para(const SplitedCscMatrix *mat, SpmvPara *para, int row_num, int col_num) {
    int chunk_num = mat->chunk_num;
    para->chunk_num = chunk_num;
    para->chunks = mat->chunks;
    CscChunk *chunk0 = mat->chunks[0];
    para->result = (double *)malloc(sizeof(double) * chunk_num * row_num);
    para->dma_over = (int *)malloc(sizeof(int) * chunk_num);
    memset(para->dma_over, 0, sizeof(int) * chunk_num);
    para->sp_col = col_num;
    para->sp_row = row_num;
    int max_block_row_num = 0;
    for (int i = 0; i < chunk0->block_num; ++i) {
        CscBlock *b = chunk0->blocks + i;
        row_num = b->row_end - b->row_begin;
        if (row_num > max_block_row_num) {
            max_block_row_num = row_num;
        }
    }
    para->max_block_row_num = max_block_row_num;
}

void spmv_para_free(SpmvPara *para) {
    free(para->dma_over);
    free(para->result);
    // TODO: free other parts
}