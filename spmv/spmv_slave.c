#include "spmv_slave.h"
#include "crts.h"
#include "slave_def.h"
#include "spmv/spmv_def.h"

__thread_local crts_rply_t get_rply = 0;
__thread_local unsigned int get_cnt = 0;

__thread_local crts_rply_t put_rply = 0;
__thread_local unsigned int put_cnt = 0;

// 双缓冲设计
typedef struct {
    void *buff[2];
    int size;
    int in_use;
} DoubleBuffering;

inline static void slave_double_buffering_new(DoubleBuffering *buff, int size) {
    buff->size = size;
    void *tmp = CRTS_pldm_malloc(size * 2);
    buff->buff[0] = tmp;
    buff->buff[1] = buff->buff[0] + size;
    buff->in_use = 1;
}

inline static void slave_double_buffering_free(DoubleBuffering *buff) {
    CRTS_pldm_free(buff->buff, buff->size * 2);
}

inline static void *slave_double_buffering_get(DoubleBuffering *buff, int size) {
    int should_use = !buff->in_use;
    buff->in_use = should_use;
    return buff->buff[should_use];
}

static void slave_csc_chunk_unpack(CscChunk *chunk) {
    int meta_seg_num = 2;
    int *meta_seg = (int *)chunk->packed_data.mem;
    int *rows_seg = (int *)chunk->packed_data.mem + meta_seg_num;
    int *col_off_seg = (int *)chunk->packed_data.mem + meta_seg[0];
    double *data_seg = (double *)chunk->packed_data.mem + meta_seg[1];

    for (int i = 0; i < chunk->block_num; ++i) {
        CscBlock *block = chunk->blocks + i;

        block->rows = rows_seg;
        rows_seg += block->data_size;

        block->data = data_seg;
        data_seg += block->data_size;

        block->col_off = col_off_seg;
        col_off_seg += block->col_num + 1;
    }
}

// 在从核上计算 csc 格式的稀疏矩阵乘向量(spmv)
void slave_csc_spmv(SpmvPara *para_mp) {
    // 这段代码中所有后缀为 _mp 的变量保存的是主存中的地址(master pointer)
    SpmvPara spmv_para;
    DMA_GET(&spmv_para, para_mp, sizeof(SpmvPara), &get_rply, get_cnt);
    int id = CRTS_tid;
    int chunk_num = spmv_para.chunk_num;
    int *dma_over_mp = spmv_para.dma_over + id;
    double *result_mp = spmv_para.result;

    // 获取所需要的 chunk 在主存中的地址
    CscChunk *chunk_mp = spmv_para.chunks[id];

    // 将所需要的 chunk 读入
    int size;
    DMA_GET(&size, &chunk_mp->size, sizeof(int), &get_rply, get_cnt);
    CscChunk *chunk = (CscChunk *)CRTS_pldm_malloc(size);
    DMA_GET(chunk, chunk_mp, size, &get_rply, get_cnt);
    int block_num = chunk->block_num;

    void *mem = CRTS_pldm_malloc(chunk->packed_data.mem_size);
    DMA_GET(mem, chunk->packed_data.mem, chunk->packed_data.mem_size, &get_rply, get_cnt);
    chunk->packed_data.mem = mem;
    slave_csc_chunk_unpack(chunk);

    CscBlock *blocks = chunk->blocks;

    // 获取该 chunk 对应 vec 的 slice
    double *vec_mp = chunk->vec;
    int col_begin = chunk->col_begin;
    int col_end = chunk->col_end;
    int slice_size = col_end - col_begin;
    double *slice = (double *)CRTS_pldm_malloc(slice_size);
    DMA_GET(slice, vec_mp + col_begin, slice_size, &get_rply, get_cnt);

    DoubleBuffering double_buff;
    slave_double_buffering_new(&double_buff, sizeof(double) * spmv_para.max_block_row_num);
    double *block_result_mp = result_mp + id;
    for (int i = 0; i < block_num; ++i) {
        CscBlock *block = blocks + i;
        int col_num = block->col_num;
        int row_num = block->row_end - block->row_begin;
        int data_size = block->data_size;
        int *rows = block->rows;
        int *col_off = block->col_off;
        double *data = block->data;

        int result_size = sizeof(double) * row_num;

        double *block_result = slave_double_buffering_get(&double_buff, result_size);

        memset(block_result, 0, sizeof(double) * row_num);

        // block 与 slice 的 spmv
        for (int j = 0; j < col_num; ++j) {
            int start = col_off[j];
            int end = col_off[j + 1];
            int num = end - start;
            double tmp = 0.;
            for (int k = 0; k < num; ++k) {
                block_result[rows[start + k]] += slice[j] * data[start + k];
            }
        }

        int last = i - 1;
        DMA_WAIT(&put_rply, put_cnt);
        DMA_IPUT(dma_over_mp, &last, sizeof(int), &put_rply, put_cnt);
        DMA_IPUT_STRIDE(block_result_mp, block_result, sizeof(double) * row_num, sizeof(double), sizeof(double) * chunk_num, &put_rply, put_cnt);

        block_result_mp += row_num * chunk_num;
    }

    int last = chunk_num - 1;
    DMA_WAIT(&put_rply, put_cnt);
    DMA_IPUT(dma_over_mp, &last, sizeof(int), &put_rply, put_cnt);

    slave_double_buffering_free(&double_buff);

    CRTS_pldm_free(slice, slice_size);
    CRTS_pldm_free(chunk->packed_data.mem, chunk->packed_data.mem_size);
    CRTS_pldm_free(chunk, size);
}