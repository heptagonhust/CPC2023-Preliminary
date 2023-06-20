#include "spmv_master.h"
#include "swperf.h"
#include <crts.h>
#include <simd.h>
#define SLAVE_CORE_NUM 64


void coo_spmv(SpmvPara *para, double *result) {
    CooChunk *chunk0 = para->chunks[0].chunk;
    int chunk_num = para->chunk_num;
    int block_num_per_chunk = chunk0->block_num;
    int *reduce_status = (int *)malloc(chunk_num * block_num_per_chunk * sizeof(int));
    memset(reduce_status, 0, chunk_num * block_num_per_chunk * sizeof(int));
    para->reduce_status = reduce_status;
    double *result_pool = para->result;
    int *dma_over = para->dma_over;
    int reduced_cnt = 0;

    CRTS_athread_spawn(slave_coo_spmv, para);
    memset(result, 0, sizeof(double) * para->sp_col);

    //! test
        unsigned long icc = 0;
        unsigned long calcu_cycles = 0;
        unsigned long retire_ins = 0;
        unsigned long calcu_ins = 0;
        int cnt = 0, cnt_1 = 0;
        penv_host0_cycle_init();
        // penv_host0_retire_inst_init();
    //! test
    while (reduced_cnt < chunk_num * block_num_per_chunk) {
        reduced_cnt = 0;
        for (int block_idx = 0; block_idx < block_num_per_chunk; ++block_idx) {
            int row_begin = chunk0->blocks[block_idx].row_begin;
            int row_num = chunk0->blocks[block_idx + 1].row_begin - chunk0->blocks[block_idx].row_begin;
            for (int chunk_idx = 0; chunk_idx < chunk_num; ++chunk_idx) {
                int flag_idx = block_num_per_chunk * chunk_idx + block_idx;
                if (dma_over[chunk_idx] >= block_idx && reduce_status[flag_idx] == 0) {

        penv_host0_cycle_count(&icc);
        unsigned long start = icc;
                    int base = chunk_num * row_begin + chunk_idx * row_num;
                    // 对界读取实现
                    int off = 0;
                    // 前部边缘处理
                    if (base % 4 != 0) {
                        off = 4 - (base % 4);
                        for (int i = 0; i < off; ++i) {
                            result[row_begin + i] += result_pool[base + i];
                        }
                    }
                    // simd
                    int round = (row_num - off) / 4;
                    int left = (row_num - off) % 4;
                    for (int i = 0; i < round; i++, off += 4) {
                        doublev4 res, pool;
                        simd_loadu(res, result + row_begin + off);
                        simd_load(pool, result_pool + base + off);
                        res = simd_vaddd(res, pool);
                        simd_storeu(res, result + row_begin + off);
                    }
                    // 后部边缘处理
                    for (int i = 0; i < left; ++i) {
                        result[row_begin + off + i] += result_pool[base + off + i];
                    }
                    reduce_status[flag_idx] = 1;
        penv_host0_cycle_count(&icc);
        calcu_cycles += icc - start;
                }
                if (reduce_status[flag_idx] == 1) reduced_cnt++;
            }
        }
    }
    //! test
        penv_host0_cycle_count(&icc);
        INFO("master total time: %d cycles, calculate time: %d cycles, calculate ins: %d\n", icc, calcu_cycles, calcu_ins);
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
    // 保证result pool 64 byte对界
    para->result = (double *)libc_aligned_malloc(sizeof(double) * chunk_num * row_num);
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
    libc_aligned_free(para->result);
}