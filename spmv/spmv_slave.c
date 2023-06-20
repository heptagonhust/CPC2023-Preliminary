#include "spmv_slave.h"
#include "crts.h"
#include "slave_def.h"
#include "swperf.h"
#include "spmv/spmv_def.h"

__thread_local crts_rply_t get_rply = 0;
__thread_local unsigned int get_cnt = 0;

__thread_local crts_rply_t put_rply = 0;
__thread_local unsigned int put_cnt = 0;


// 双缓冲设计
// typedef struct {
//     void *buff[2];
//     int size;
//     int in_use;
// } DoubleBuffering;

// inline static void slave_double_buffering_new(DoubleBuffering *buff, int size) {
//     buff->size = size;
//     void *tmp = CRTS_pldm_malloc(size * 2);
//     buff->buff[0] = tmp;
//     buff->buff[1] = buff->buff[0] + size;
//     buff->in_use = 1;
// }

// inline static void *slave_double_buffering_get(DoubleBuffering *buff) {
//     int should_use = !buff->in_use;
//     buff->in_use = should_use;
//     return buff->buff[should_use];
// }

// inline static void slave_double_buffering_free(DoubleBuffering *buff) {
//     CRTS_pldm_free(buff->buff[0], buff->size * 2);
// }

#define MAX_CELL_NUM 88362
// __thread_local_share double z[MAX_CELL_NUM] __attribute__((aligned(64)));
// __thread_local_share double p[MAX_CELL_NUM] __attribute__((aligned(64)));
__thread_local_share double vec[MAX_CELL_NUM] __attribute__((aligned(64)));


void slave_coo_spmv(SpmvPara *para_mp) {
    // 这段代码中所有后缀为 _mp 的变量保存的是主存中的地址(master pointer)
    SpmvPara spmv_para;
    DMA_GET(&spmv_para, para_mp, sizeof(SpmvPara), &get_rply, get_cnt);
    int id = CRTS_tid;
    int chunk_num = spmv_para.chunk_num;

    // 拉取chunk元数据
    CooChunk *chunk = (CooChunk *)CRTS_pldm_malloc(spmv_para.chunks[id].mem_size);
    DMA_GET(chunk, spmv_para.chunks[id].chunk, spmv_para.chunks[id].mem_size, &get_rply, get_cnt);

    // Vector space
    int rows = chunk->row_end - chunk->row_begin;
    CRTS_memcpy_sldm(&vec, spmv_para.vec, spmv_para.sp_col * sizeof(double), MEM_TO_LDM);

    double *result_ldm = (double *)CRTS_pldm_malloc(rows * sizeof(double));
    memset(result_ldm, 0, rows * sizeof(double));
    int chunk_data_size = chunk->blocks[chunk->block_num].block_off;
    if (chunk_data_size % 2) chunk_data_size++;
    int ldm_left_size = CRTS_pldm_get_free_size();

    // ! just for debug
    // 单次传入即可完成
    if (ldm_left_size >= chunk_data_size * sizeof(double) * 1.5 + 64 * 3) {
        double *data = (double *)CRTS_pldm_malloc(chunk_data_size * sizeof(double));
        DMA_IGET(data, chunk->data, chunk_data_size * sizeof(double), &get_rply, get_cnt);
        uint16_t *row_idx = (uint16_t *)CRTS_pldm_malloc(chunk_data_size * sizeof(uint16_t));
        DMA_IGET(row_idx, chunk->row_idx, chunk_data_size * sizeof(uint16_t), &get_rply, get_cnt);
        uint16_t *col_idx = (uint16_t *)CRTS_pldm_malloc(chunk_data_size * sizeof(uint16_t));
        DMA_IGET(col_idx, chunk->col_idx, chunk_data_size * sizeof(uint16_t), &get_rply, get_cnt);
        DMA_WAIT(&get_rply, get_cnt);

        // block processing loop with double buffering
        for (int i = 0; i < chunk->block_num; ++i) {
            int block_data_size = chunk->blocks[i+1].block_off - chunk->blocks[i].block_off;
            if (block_data_size != 0) {
                int block_data_offset = chunk->blocks[i].block_off;
                int block_col_begin = chunk->blocks[i].col_begin;
                for (int j = 0; j < block_data_size; ++j) {
                    result_ldm[row_idx[block_data_offset + j]] += data[block_data_offset + j] * vec[block_col_begin + col_idx[block_data_offset + j]];
                }
            }
        }
        DMA_IPUT(spmv_para.result + chunk->row_begin, result_ldm, rows * sizeof(double), &put_rply, put_cnt);
        CRTS_pldm_free(data, chunk_data_size * sizeof(double));
        CRTS_pldm_free(row_idx, chunk_data_size * sizeof(uint16_t));
        CRTS_pldm_free(col_idx, chunk_data_size * sizeof(uint16_t));
    }

    else {
        printf("[ERROR] Coo matrix size is more than expected!\n");
    }

    CRTS_pldm_free(chunk, spmv_para.chunks[id].mem_size);
    DMA_WAIT(&put_rply, put_cnt);
    CRTS_pldm_free(result_ldm, rows * sizeof(double));
}