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
    CRTS_pldm_free(buff->buff[0], buff->size * 2);
}


void slave_coo_spmv(SpmvPara *para_mp) {
    // profile
    unsigned long icc = 0;
    unsigned long wait_cycle = 0;
    unsigned long signal_cycle = 0;
    penv_slave0_cycle_init();

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
        if (id == 0) {
printf("LINE 80\n");
        }
        int last_idma_round = -1;
        int *dma_over_mem = spmv_para.dma_over + id;
        for (int i = 0; i < chunk->block_num; ++i) {
            int row_num = chunk->blocks[i+1].row_begin - chunk->blocks[i].row_begin;
            int col_num = chunk->blocks[i].col_num;
            int block_data_offset = chunk->blocks[i].block_off;
            int block_data_size = chunk->blocks[i+1].block_off - chunk->blocks[i].block_off;
            if (block_data_size == 0) {
        if (id == 0) {
printf("LINE 89\n");
        }
                int *reduce_flag = spmv_para.reduce_status + i * chunk_num + id;
                    //! test
                    penv_slave0_cycle_count(&icc);
                    unsigned long start_cycle = icc;
                    //! test
                CRTS_ssig_host(reduce_flag);
                int flag = 1;
                DMA_IPUT(reduce_flag, &flag, sizeof(int), &put_rply, put_cnt);
                    //! test
                    penv_slave0_cycle_count(&icc);
                    unsigned long end_cycle = icc;
                    signal_cycle += end_cycle - start_cycle;
                    //! test
        if (id == 0) {
printf("LINE 92\n");
        }
            }
            else {
        if (id == 0) {
printf("LINE 95\n");
        }
                double *block_result_ldm = slave_double_buffering_get(&double_buff);
                memset(block_result_ldm, 0, row_num * sizeof(double));
                double *block_result_mem = spmv_para.result + spmv_para.chunk_num * chunk->blocks[i].row_begin + row_num * id;

                for (int j = 0; j < block_data_size; ++j) {
                    block_result_ldm[row_idx[block_data_offset + j]] += data[block_data_offset + j] * vec[col_idx[block_data_offset + j]];
                }
            
                last_idma_round = i - 1;
                    penv_slave0_cycle_count(&icc);
                    unsigned long start_cycle = icc;
                DMA_WAIT(&put_rply, put_cnt);
                    penv_slave0_cycle_count(&icc);
                    unsigned long end_cycle = icc;
                    wait_cycle += end_cycle - start_cycle;
                DMA_IPUT(dma_over_mem, &last_idma_round, sizeof(int), &put_rply, put_cnt);
                DMA_IPUT(block_result_mem, block_result_ldm, sizeof(double) * row_num, &put_rply, put_cnt);
        if (id == 0) {
printf("LINE 113\n");
        }
            }
        }
        last_idma_round = chunk_num - 1;
            penv_slave0_cycle_count(&icc);
            unsigned long start_cycle = icc;
        DMA_WAIT(&put_rply, put_cnt);
            penv_slave0_cycle_count(&icc);
            unsigned long end_cycle = icc;
            wait_cycle += end_cycle - start_cycle;
        DMA_IPUT(dma_over_mem, &last_idma_round, sizeof(int), &put_rply, put_cnt);
        CRTS_pldm_free(data, chunk_data_size * sizeof(double));
        CRTS_pldm_free(row_idx, chunk_data_size * sizeof(uint16_t));
        CRTS_pldm_free(col_idx, chunk_data_size * sizeof(uint16_t));
        DMA_WAIT(&put_rply, put_cnt);
        if (id == 0) {
printf("LINE 128\n");
        }
    }

    else {
        printf("[ERROR] Coo matrix size is more than expected!\n");
    }

    slave_double_buffering_free(&double_buff);
    CRTS_pldm_free(vec, cols * sizeof(double));
    CRTS_pldm_free(chunk, spmv_para.chunks[id].mem_size);
        if (id == 0) {
printf("LINE 138\n");
        }
    // ldm_left_size = CRTS_pldm_get_free_size();

    penv_slave0_cycle_count(&icc);
    printf("Slave core %d: Total: %lld cycles, DMA wait: %lld cycles, signal : %lld cycles\n", id, icc, wait_cycle, signal_cycle);
}