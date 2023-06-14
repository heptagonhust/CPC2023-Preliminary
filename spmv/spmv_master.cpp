#include "spmv_master.h"

#define SLAVE_CORE_NUM 64
void csc_spmv(SpmvPara *para, double *result) {
    CRTS_athread_spawn(slave_csc_spmv, &para);

    memset(result, 0, sizeof(double) * para->sp_row);

    int chunk_num = para->chunk_num;
    double *result_pool = para->result;
    int *dma_over = para->dma_over;
    CscChunk *chunk0 = para->chunks[0];
    int block_num_per_chunk = chunk0->block_num;

    bool *chunk_over = (bool *)malloc(sizeof(bool) * chunk_num);
    memset(chunk_over, 0, sizeof(bool) * 64);
    for (int i = 0; i < block_num_per_chunk; ++i) {
        CscBlock *block = chunk0->blocks + i;
        int row_begin = block->row_begin;
        int row_num = block->row_end - row_begin;

        double *result_slice = result + row_begin;
        double *chunk_result_base = result_pool + chunk_num * row_begin;
        for (;;) {
            bool get_out = true;
            for (int j = 0; j < chunk_num; ++j) {
                if (dma_over[j] >= i) {
                    double *chunk_result = chunk_result_base + j;
                    for (int k = 0; k < row_num; ++k) {
                        result_slice[k] += *(chunk_result + k * chunk_num);
                    }
                    chunk_over[j] = true;
                }
                get_out &= chunk_over[j];
            }

            if (get_out) {
                break;
            }
        }

        memset(chunk_over, 0, sizeof(bool) * 64);
    }

    free(chunk_over);
}

void splited_csc_matrix_to_spmv_para(
    const SplitedCscMatrix *mat,
    SpmvPara *para,
    int row_num,
    int col_num) {
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