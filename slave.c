#include <slave.h>
#include <crts.h>
#include <math.h>
#include <simd.h>
#include "vector_def.h"
#include "spmv_def.h"
#include "spmv_slave.h"
#include "slave_def.h"
#include "perf.h"

__thread_local crts_rply_t get_rply = 0;
__thread_local unsigned int get_cnt = 0;

__thread_local crts_rply_t put_rply = 0;
__thread_local unsigned int put_cnt = 0;

// #define MAX_CELL_NUM 88362
#define REDUCE_BUF_SIZE 64
#define BLOCK_SIZE 64
// __thread_local_share double p[MAX_CELL_NUM] __attribute__((aligned(64)));
// __thread_local_share double z[MAX_CELL_NUM] __attribute__((aligned(64)));
__thread_local uint16_t non_0_block_num = 0;
__thread_local double reduce_buf[REDUCE_BUF_SIZE] __attribute__((aligned(64)));
__thread_local uint16_t non_0_block_idx[BLOCK_SIZE] __attribute__((aligned(64)));
__thread_local_share double* z_list[BLOCK_SIZE] __attribute__((aligned(64)));
__thread_local_share double* p_list[BLOCK_SIZE] __attribute__((aligned(64)));

__thread_local PerfEnv env;


inline void get_non_0_blocks(CooChunk *chunk) {
    non_0_block_num = 0;
    int id = CRTS_tid;
    for (int i = id; i < chunk->block_num; ++i) {
        if (chunk->blocks[i + 1].block_off - chunk->blocks[i].block_off != 0) {
            non_0_block_idx[non_0_block_num] = i;
            ++non_0_block_num;
        }
    }
    for (int i = 0; i < id; ++i) {
        if (chunk->blocks[i + 1].block_off - chunk->blocks[i].block_off != 0) {
            non_0_block_idx[non_0_block_num] = i;
            ++non_0_block_num;
        }
    }
}

//TODO can be optimized
inline static void slave_precondition_coo(CooChunk *chunk, double *M, double *M_1, double *r, double *z, double *g, int vec_num, int vec_begin, DoubleBuffering *buff) {
    slave_Mul(z, M_1, r, vec_num);
    for (int deg = 1; deg < 2; deg++) {
        CRTS_ssync_array();
slave_perf_begin(&env, PERF_INCLUSIVE, "SPMV");
        slave_coo_spmv(chunk, z, z_list, g, buff, non_0_block_num, &non_0_block_idx[0]);
slave_perf_end(&env, "SPMV");
        slave_MulSub(g, z, M, vec_num);
        slave_SubMul(z, r, g, M_1, vec_num);
    }
}

inline static void slave_init_M_1(double *M, double *M_1, int vec_num) {
    int i;
    doublev8 M8, M_18;
    for (i = 0; i < vec_num - 8; i += 8) {
        simd_load(M8, M + i);
        M_18 = simd_vfrecpd(M8);
        simd_store(M_18, M_1 + i);
    }

    for (; i < vec_num; ++i) {
        M_1[i] = 1.0 / M[i];
    }
}

void slave_MainLoop(MainLoopPara *para_mem) {
slave_init_perf_env(&env, 0);
slave_perf_begin(&env, PERF_INCLUSIVE, "Main Loop");
slave_perf_begin(&env, PERF_EXCLUSIVE, "DMA GET DATA");
    int id = CRTS_tid;
    
    // get meta data
    MainLoopPara para;
    DMA_GET(&para, para_mem, sizeof(MainLoopPara), &get_rply, get_cnt);
    int vec_begin = para.ntask[id].col_start;
    int vec_num = para.ntask[id].col_num;
    SpmvPara spmv_para;
    DMA_GET(&spmv_para, para.spmv_para, sizeof(SpmvPara), &get_rply, get_cnt);
    double *M = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *M_1 = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    DMA_IGET(M, para.M + vec_begin, vec_num * sizeof(double), &get_rply, get_cnt);

    CooChunk *chunk = (CooChunk *)CRTS_pldm_malloc(spmv_para.chunks[id].mem_size);
    DMA_IGET(chunk, spmv_para.chunks[id].chunk, spmv_para.chunks[id].mem_size, &get_rply, get_cnt);
    DMA_WAIT(&get_rply, get_cnt);

    DoubleBuffering buffer;
    slave_double_buffering_new(&buffer, spmv_para.max_block_row_num);

    //* get vector data
    double *r = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *g = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *x = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *Ax = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *z = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    z_list[id] = z;
    double *p = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    p_list[id] = p;
    DMA_IGET(r, para.r + vec_begin, vec_num * sizeof(double), &get_rply, get_cnt);
    DMA_IGET(x, para.x + vec_begin, vec_num * sizeof(double), &get_rply, get_cnt);
    DMA_IGET(p, para.p + vec_begin, vec_num * sizeof(double), &get_rply, get_cnt);

    get_non_0_blocks(chunk);
    slave_init_M_1(M, M_1, vec_num);

    DMA_WAIT(&get_rply, get_cnt);

    //* get spmv chunk data
    int chunk_data_size = chunk->blocks[chunk->block_num].block_off;
    if (chunk_data_size % 2) chunk_data_size++;
    //! only for debug
    int ldm_left_size = CRTS_pldm_get_free_size();
    if (ldm_left_size >= chunk_data_size * sizeof(double) * 1.5 + 64 * 3) {
        double *data = (double *)CRTS_pldm_malloc(chunk_data_size * sizeof(double));
        DMA_IGET(data, chunk->data, chunk_data_size * sizeof(double), &get_rply, get_cnt);
        uint16_t *row_idx = (uint16_t *)CRTS_pldm_malloc(chunk_data_size * sizeof(uint16_t));
        DMA_IGET(row_idx, chunk->row_idx, chunk_data_size * sizeof(uint16_t), &get_rply, get_cnt);
        uint16_t *col_idx = (uint16_t *)CRTS_pldm_malloc(chunk_data_size * sizeof(uint16_t));
        DMA_IGET(col_idx, chunk->col_idx, chunk_data_size * sizeof(uint16_t), &get_rply, get_cnt);
        DMA_WAIT(&get_rply, get_cnt);
        chunk->data = data;
        chunk->row_idx = row_idx;
        chunk->col_idx = col_idx;
    }
    else {
        printf("[ERROR] Coo matrix size is more than expected!\n");
    }
slave_perf_end(&env, "DMA GET DATA");

slave_perf_begin(&env, PERF_INCLUSIVE, "DATA PROCESSING");
    CRTS_ssync_array();
slave_perf_begin(&env, PERF_INCLUSIVE, "SPMV");
    slave_coo_spmv(chunk, p, p_list, Ax, &buffer, non_0_block_num, &non_0_block_idx[0]);
slave_perf_end(&env, "SPMV");
    for (int i = 0; i < vec_num; ++i) {
        r[i] = r[i] - Ax[i];
    }
    double alpha, beta;
    double residual = slave_Reduce(r, vec_num, &reduce_buf[0]);
    double init_residual = residual;
    double sumprod, sumprod_old;

    int iter = 0;
    if (fabs(residual / para.normfactor) > para.tolerance) {
        do {
            if (iter == 0) {
                slave_precondition_coo(chunk, M, M_1, r, z, g, vec_num, vec_begin, &buffer);
                sumprod = slave_MulReduceZR(r, z, vec_num, &reduce_buf[0]);
                memcpy(p, z, vec_num * sizeof(double));
            }
            else {
                sumprod_old = sumprod;
                slave_precondition_coo(chunk, M, M_1, r, z, g, vec_num, vec_begin, &buffer);
                sumprod = slave_MulReduceZR(r, z, vec_num, &reduce_buf[0]);
                beta = sumprod / sumprod_old;
                slave_MulAdd(p, z, beta, vec_num, vec_begin);
            }
            CRTS_ssync_array();
slave_perf_begin(&env, PERF_INCLUSIVE, "SPMV");
            slave_coo_spmv(chunk, p, p_list, Ax, &buffer, non_0_block_num, &non_0_block_idx[0]);
slave_perf_end(&env, "SPMV");
            alpha = sumprod / slave_MulReducepAx(p, Ax, vec_num, &reduce_buf[0]);
            slave_Updatexr(x, r, p, Ax, alpha, vec_num);
            residual = slave_Reduce(r, vec_num, &reduce_buf[0]);
        } while (++iter < para.maxIter && (residual / para.normfactor) >= para.tolerance);
    }
    slave_perf_end(&env, "DATA PROCESSING");
    DMA_IPUT(para.init_residual, &init_residual, sizeof(double), &put_rply, put_cnt);
    DMA_IPUT(para.final_residual, &residual, sizeof(double), &put_rply, put_cnt);
    DMA_IPUT(para.iter, &iter, sizeof(int), &put_rply, put_cnt);

    slave_double_buffering_free(&buffer);
    CRTS_pldm_free(chunk->data, chunk_data_size * sizeof(double));
    CRTS_pldm_free(chunk->row_idx, chunk_data_size * sizeof(uint16_t));
    CRTS_pldm_free(chunk->col_idx, chunk_data_size * sizeof(uint16_t));
    CRTS_pldm_free(r, vec_num * sizeof(double));
    CRTS_pldm_free(g, vec_num * sizeof(double));
    CRTS_pldm_free(x, vec_num * sizeof(double));
    CRTS_pldm_free(p, vec_num * sizeof(double));
    CRTS_pldm_free(z, vec_num * sizeof(double));
    CRTS_pldm_free(Ax, vec_num * sizeof(double));
    CRTS_pldm_free(M, vec_num * sizeof(double));
    CRTS_pldm_free(M_1, vec_num * sizeof(double));
    CRTS_pldm_free(chunk, spmv_para.chunks[id].mem_size);
    DMA_WAIT(&put_rply, put_cnt);
slave_perf_end(&env, "Main Loop");
if (id == 0) slave_perf_report(&env);
}