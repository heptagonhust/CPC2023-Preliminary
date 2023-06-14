#include "spmv_slave.h"
#include "crts.h"
#include "slave_def.h"
#include "spmv/spmv_def.h"

__thread_local crts_rply_t get_rply = 0;
__thread_local unsigned int get_cnt = 0;

__thread_local crts_rply_t put_rply = 0;
__thread_local unsigned int put_cnt = 0;

#define BUF_ITEM 1024
__thread_local double redu_buf[BUF_ITEM] __attribute__((aligned(64)));

void slave_double_buffering_free(DoubleBuffering *buff) {
    for (int i = 0; i < 2; ++i) {
        void *b = buff->buff[i];
        int s = buff->size[i];
        if (s) {
            CRTS_pldm_free(b, s);
        }
    }
}

void *slave_double_buffering_get(DoubleBuffering *buff, int size) {
    int should_use = !buff->in_use;
    if (buff->size[should_use] < size) {
        CRTS_pldm_free(buff->buff[should_use], buff->size[should_use]);
        buff->buff[should_use] = CRTS_pldm_malloc(size);
        buff->size[should_use] = size;
    }
    buff->in_use = should_use;
    return buff->buff[should_use];
}

void slave_csc_spmv(void *para) {
    SpmvPara spmv_para;
    DMA_GET(&spmv_para, para, sizeof(SpmvPara), &get_rply, get_cnt);
    int id = CRTS_tid;
    int chunk_num = spmv_para.chunk_num;
    int *dma_over = spmv_para.dma_over + id;
    double *result = spmv_para.result;

    int col = spmv_para.sp_col;
    int row = spmv_para.sp_row;

    // 获取所需要的 chunk 在主存中的地址
    CscChunk **chunk_ptr = spmv_para.chunks + id;

    // 将所需要的 chunk 读入
    int size;
    DMA_GET(&size, &((CscChunk *)chunk_ptr)->size, sizeof(int), &get_rply, get_cnt);
    CscChunk *chunk = (CscChunk *)CRTS_pldm_malloc(size);
    DMA_GET(chunk, chunk_ptr, size, &get_rply, get_cnt);
    int block_num = chunk->block_num;

    // 一次性读入所有的 blocks
    CscBlock *blocks = (CscBlock *)CRTS_pldm_malloc(sizeof(CscBlock) * block_num);
    for (int i = 0; i < block_num; ++i) {
        CscBlock *block_ptr = chunk->blocks[i];
        CscBlock *b = blocks + i;
        DMA_GET(b, block_ptr, sizeof(CscBlock), &get_rply, get_cnt);
        int data_size = b->data_size;
        int total_size = b->total_size;
        double *buffer = (double *)CRTS_pldm_malloc(total_size);
        b->data = (double *)buffer;
        b->rows = (int *)(b->data + data_size);
        b->col_off = b->rows + data_size;
    }

    // 获取该 chunk 对应 vec 的 slice
    double *vec = chunk->vec;
    int col_begin = chunk->col_begin;
    int col_end = chunk->col_end;
    int slice_size = col_end - col_begin;
    double *slice = (double *)CRTS_pldm_malloc(slice_size);
    DMA_GET(slice, vec + col_begin, slice_size, &get_rply, get_cnt);

    DoubleBuffering double_buff;
    double *block_result_ptr = result + id;
    for (int i = 0; i < block_num; ++i) {
        CscBlock *block = blocks + i;
        int col_num = block->col_num;
        int data_size = block->data_size;
        int *rows = block->rows;
        int *col_off = block->col_off;
        double *data = block->data;

        int result_size = sizeof(double) * col_num;

        // 双缓冲并每次比较一下，减少分配次数
        double *block_result = slave_double_buffering_get(&double_buff, result_size);

        // block 与 slice 的 spmv
        for (int j = 0; j < col_num; ++j) {
            int start = col_off[j];
            int end = col_off[j + 1];
            int num = end - start;
            double tmp = 0.;
            for (int k = 0; k < num; ++k) {
                tmp += slice[rows[start + k]] * data[start + k];
            }
            block_result[j] = tmp;
        }

        int last = i - 1;
        DMA_IPUT(dma_over, &last, sizeof(int), &put_rply, put_cnt);
        DMA_WAIT(&put_rply, put_cnt);
        DMA_IPUT_STRIDE(block_result_ptr, block_result, sizeof(double) * col_num, sizeof(double), sizeof(double) * chunk_num, &put_rply, put_cnt);

        block_result_ptr += col_num * chunk_num;
    }

    int last = chunk_num - 1;
    DMA_IPUT(dma_over, &last, sizeof(int), &put_rply, put_cnt);

    slave_double_buffering_free(&double_buff);

    for (int i = 0; i < block_num; ++i) {
        int total_size = blocks->total_size;
        CRTS_pldm_free(blocks + i, total_size);
    }
    CRTS_pldm_free(blocks, sizeof(CscBlock *) * block_num);

    CRTS_pldm_free(slice, slice_size);
    CRTS_pldm_free(chunk, size);
}