#include "spmv_master.h"
#include "swperf.h"
#include <crts.h>

#define SLAVE_CORE_NUM 64


void coo_spmv(SpmvPara *para, double *result) {
    CooChunk *chunk0 = para->chunks[0].chunk;
    int chunk_num = para->chunk_num;
    int block_num_per_chunk = chunk0->block_num;
    int *reduce_status = (int *)calloc(chunk_num * block_num_per_chunk, sizeof(int));
    para->reduce_status = reduce_status;
    double *result_pool = para->result;
    int *dma_over = para->dma_over;
    int reduced_cnt = 0;

    CRTS_athread_spawn(slave_coo_spmv, para);
    memset(result, 0, sizeof(double) * para->sp_col);

    //! test
        unsigned long icc;
        penv_host0_cycle_init();
    //! test
    while (reduced_cnt < chunk_num * block_num_per_chunk) {
        reduced_cnt = 0;
        for (int block_idx = 0; block_idx < block_num_per_chunk; ++block_idx) {
            int row_begin = chunk0->blocks[block_idx].row_begin;
            int row_num = chunk0->blocks[block_idx + 1].row_begin - chunk0->blocks[block_idx].row_begin;
            for (int chunk_idx = 0; chunk_idx < chunk_num; ++chunk_idx) {
                int flag_idx = chunk_num * block_idx + chunk_idx;
                if (dma_over[chunk_idx] >= block_idx && !reduce_status[flag_idx]) {
                    for (int off = 0; off < row_num; ++off) {
                        result[row_begin + off] += result_pool[chunk_num * row_begin + chunk_idx * row_num + off];
                    }
                    reduce_status[flag_idx] = 1;
                }
                if (reduce_status[flag_idx] == 1) ++reduced_cnt;
            }
        }
    }
    //! test
        penv_host0_cycle_count(&icc);
        INFO("master reduce time: %d cycles\n", icc);
    //! test

    free(reduce_status);
    CRTS_athread_join();
}

void spmv_para_from_splited_coo_matrix(
    const SplitedCooMatrix *mat,
    SpmvPara *para,
    double *vec,
    int row_num,
    int col_num) {
    int chunk_num = mat->chunk_num;
    para->chunk_num = chunk_num;
    para->sp_row = row_num;
    para->sp_col = col_num;
    para->result = (double *)malloc(sizeof(double) * chunk_num * row_num);
    memset(para->result, 0, sizeof(double) * chunk_num * row_num);
    para->dma_over = (int *)malloc(sizeof(int) * chunk_num);
    memset(para->dma_over, 0, sizeof(int) * chunk_num);
    para->vec = vec;
    memcpy(para->chunks, mat->chunks, SLAVE_CORE_NUM * sizeof(SizedCooChunk));
    // all chunks have the same block size
    CooChunk *chunk0 = mat->chunks[0].chunk;
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