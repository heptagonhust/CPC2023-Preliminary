#include <slave.h>
#include <crts.h>
#include <math.h>
#include "vector_def.h"
#include "spmv_def.h"
#include "spmv_slave.h"
#include "slave_def.h"

__thread_local crts_rply_t get_rply = 0;
__thread_local unsigned int get_cnt = 0;

__thread_local crts_rply_t put_rply = 0;
__thread_local unsigned int put_cnt = 0;

#define MAX_CELL_NUM 88362
#define REDUCE_BUF_SIZE 64
__thread_local_share double p[MAX_CELL_NUM] __attribute__((aligned(64)));
__thread_local_share double z[MAX_CELL_NUM] __attribute__((aligned(64)));
__thread_local double reduce_buf[REDUCE_BUF_SIZE] __attribute__((aligned(64)));
__thread_local int cells;

//TODO can be optimized
void slave_precondition_coo(CooChunk *chunk, double *M, double *M_1, double *r, double *z, double *g, int vec_num, int vec_begin) {
    slave_Mul(z + vec_begin, M_1, r, vec_num);
    for (int deg = 1; deg < 2; deg++) {
        CRTS_ssync_sldm();
        slave_coo_spmv(chunk, z, g);
        slave_MulSub(g, z + vec_begin, M, vec_num);
        CRTS_ssync_sldm();
        slave_SubMul(z + vec_begin, r, g, M_1, vec_num);
    }
}

void slave_MainLoop(MainLoopPara *para_mem) {
    int id = CRTS_tid;
    //* get meta data
    MainLoopPara para;
    DMA_GET(&para, para_mem, sizeof(MainLoopPara), &get_rply, get_cnt);
    int vec_begin = para.ntask[id].col_start;
    int vec_num = para.ntask[id].col_num;
    SpmvPara spmv_para;
    DMA_GET(&spmv_para, para.spmv_para, sizeof(SpmvPara), &get_rply, get_cnt);

    cells = spmv_para.sp_col;

    CooChunk *chunk = (CooChunk *)CRTS_pldm_malloc(spmv_para.chunks[id].mem_size);
    DMA_IGET(chunk, spmv_para.chunks[id].chunk, spmv_para.chunks[id].mem_size, &get_rply, get_cnt);
    //* get vector data
    double *r = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *g = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *x = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *Ax = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *M = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    double *M_1 = (double *)CRTS_pldm_malloc(vec_num * sizeof(double));
    DMA_IGET(r, para.r + vec_begin, vec_num * sizeof(double), &get_rply, get_cnt);
    DMA_IGET(x, para.x + vec_begin, vec_num * sizeof(double), &get_rply, get_cnt);
    DMA_IGET(M, para.M + vec_begin, vec_num * sizeof(double), &get_rply, get_cnt);
    DMA_IGET(M_1, para.M_1 + vec_begin, vec_num * sizeof(double), &get_rply, get_cnt);
    CRTS_memcpy_sldm(&p, para.p, spmv_para.sp_col * sizeof(double), MEM_TO_LDM);
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

    slave_coo_spmv(chunk, p, Ax);
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
                slave_precondition_coo(chunk, M, M_1, r, z, g, vec_num, vec_begin);
                sumprod = slave_MulReduceZR(r, z + vec_begin, vec_num, &reduce_buf[0]);
                memcpy(p, z, spmv_para.sp_col * sizeof(double));
            }
            else {
                sumprod_old = sumprod;
                slave_precondition_coo(chunk, M, M_1, r, z, g, vec_num, vec_begin);
                sumprod = slave_MulReduceZR(r, z + vec_begin, vec_num, &reduce_buf[0]);
                beta = sumprod / sumprod_old;
                slave_MulAdd(p + vec_begin, z + vec_begin, beta, vec_num);
            }
            CRTS_ssync_array();
            slave_coo_spmv(chunk, p, Ax);
            alpha = sumprod / slave_MulReducepAx(p + vec_begin, Ax, vec_num, &reduce_buf[0]);
            slave_Updatexr(x, r, p + vec_begin, Ax, alpha, vec_num);
            residual = slave_Reduce(r, vec_num, &reduce_buf[0]);
        } while (++iter < para.maxIter && (residual / para.normfactor) >= para.tolerance);
    }
    DMA_IPUT(para.init_residual, &init_residual, sizeof(double), &put_rply, put_cnt);
    DMA_IPUT(para.final_residual, &residual, sizeof(double), &put_rply, put_cnt);
    DMA_IPUT(para.iter, &iter, sizeof(int), &put_rply, put_cnt);

    CRTS_pldm_free(chunk->data, chunk_data_size * sizeof(double));
    CRTS_pldm_free(chunk->row_idx, chunk_data_size * sizeof(uint16_t));
    CRTS_pldm_free(chunk->col_idx, chunk_data_size * sizeof(uint16_t));
    CRTS_pldm_free(r, vec_num * sizeof(double));
    CRTS_pldm_free(g, vec_num * sizeof(double));
    CRTS_pldm_free(x, vec_num * sizeof(double));
    CRTS_pldm_free(Ax, vec_num * sizeof(double));
    CRTS_pldm_free(M, vec_num * sizeof(double));
    CRTS_pldm_free(M_1, vec_num * sizeof(double));
    CRTS_pldm_free(chunk, spmv_para.chunks[id].mem_size);
    DMA_WAIT(&put_rply, put_cnt);
}