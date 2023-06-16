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

inline static void *slave_double_buffering_get(DoubleBuffering *buff) {
    int should_use = !buff->in_use;
    buff->in_use = should_use;
    return buff->buff[should_use];
}

inline static void slave_double_buffering_free(DoubleBuffering *buff) {
    CRTS_pldm_free(buff->buff, buff->size * 2);
}


void slave_csc_spmv(SpmvPara *para_mp) {
    // 这段代码中所有后缀为 _mp 的变量保存的是主存中的地址(master pointer)
    SpmvPara spmv_para;
    DMA_GET(&spmv_para, para_mp, sizeof(SpmvPara), &get_rply, get_cnt);
    int id = CRTS_tid;
    int chunk_num = spmv_para.chunk_num;

    // 拉取chunk元数据
    CooChunk *chunk = (CooChunk *)CRTS_pldm_malloc(spmv_para.chunks[id].mem_size);
    DMA_GET(chunk, spmv_para.chunks[id].chunk, spmv_para.chunks[id].mem_size, &get_rply, get_cnt);

    // 元数据DMA完成，进行内存容量分析
    // 先分配必须的空间
    // Double Buffer space
    DoubleBuffering double_buff;
    slave_double_buffering_new(&double_buff, sizeof(double) * spmv_para.max_block_row_num);
    // Vector space
    int cols = chunk->col_end - chunk->col_begin;
    double *vec = (double *)CRTS_pldm_malloc(cols * sizeof(double)); 
    DMA_GET(vec, spmv_para.vec + chunk->col_begin, cols * sizeof(double), &get_rply, get_cnt);
    int chunk_data_size = chunk->blocks[chunk->block_num].block_off - chunk->blocks[0].block_off;
    int ldm_left_size = CRTS_pldm_get_free_size();

    // chunk data
    double *data;
    int *row_idx, *col_idx;
    // 单次传入即可完成
    if (ldm_left_size >= chunk_data_size * sizeof(double) * 2 + 64 * 3) {
        data = (double *)CRTS_pldm_malloc(chunk_data_size * sizeof(double));
        DMA_IGET(data, chunk->data, chunk_data_size * sizeof(double), &get_rply, get_cnt);
        row_idx = (int *)CRTS_pldm_malloc(chunk_data_size * sizeof(int));
        DMA_IGET(row_idx, chunk->row_idx, chunk_data_size * sizeof(int), &get_rply, get_cnt);
        col_idx = (int *)CRTS_pldm_malloc(chunk_data_size * sizeof(int));
        DMA_IGET(col_idx, chunk->col_idx, chunk_data_size * sizeof(int), &get_rply, get_cnt);
        DMA_WAIT(&get_rply, get_cnt);

        // block processing loop with double buffering
        int last_idma_round = -1;
        int *dma_over_mem = spmv_para.dma_over + id;
        double *block_result_mem = spmv_para.result + id;
        for (int i = 0; i < chunk->block_num; ++i) {
            int row_num = chunk->blocks[i+1].row_begin - chunk->blocks[i].row_begin;
            int col_num = chunk->blocks[i].col_num;
            int block_data_size = chunk->blocks[i+1].block_off - chunk->blocks[i].block_off;
            double *block_result_ldm = slave_double_buffering_get(&double_buff);
            memset(block_result_ldm, 0, row_num * sizeof(double));

            for (int j = 0; j < block_data_size; ++j) {
                block_result_ldm[row_idx[i]] += data[i] * vec[col_idx[i]];
            }
            
            last_idma_round = i - 1;
            DMA_WAIT(&put_rply, put_cnt);
            DMA_IPUT(dma_over_mp, &last, sizeof(int), &put_rply, put_cnt);
            DMA_IPUT_STRIDE(block_result_mp, block_result, sizeof(double) * row_num, sizeof(double), sizeof(double) * chunk_num, &put_rply, put_cnt);

            block_result_mem += row_num * chunk_num;
        }
        last_idma_round = chunk_num - 1;
        DMA_WAIT(&put_rply, put_cnt);
        DMA_IPUT(dma_over_mp, &last, sizeof(int), &put_rply, put_cnt);
        
    }
    // 需通过多次iDMA传输完成
    else {

    }


    slave_double_buffering_free(&double_buff);
    CRTS_pldm_free(vec, cols * sizeof(double));
    CRTS_pldm_free(data, chunk_data_size * sizeof(double));
    CRTS_pldm_free(row_idx, chunk_data_size * sizeof(int));
    CRTS_pldm_free(col_idx, chunk_data_size * sizeof(int));
    DMA_WAIT(&put_rply, put_cnt);
}