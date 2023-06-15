#include "vector_slave.h"

#include <crts.h>
#include <slave.h>

#include "pcg_def.h"
#include "vector_def.h"

#define reduceBufferSize 64
__thread_local crts_rply_t DMARply = 0;
__thread_local unsigned int DMARplyCount = 0;
__thread_local double reducebuf[reduceBufferSize] __attribute__((aligned(64)));
// #define dataBufferSize 1408
// __thread_local double p[dataBufferSize] __attribute__((aligned(64)));
// __thread_local double z[dataBufferSize] __attribute__((aligned(64)));
// __thread_local double M[dataBufferSize] __attribute__((aligned(64)));
// __thread_local double r[dataBufferSize] __attribute__((aligned(64)));
// __thread_local double g[dataBufferSize] __attribute__((aligned(64)));
// __thread_local double Ax[dataBufferSize] __attribute__((aligned(64)));
// __thread_local double x[dataBufferSize] __attribute__((aligned(64)));

// --------------------------------------------------------------------------------

void slave_MulAdd(MulAddPara *para) {
    MulAddPara slavePara;
    //接收结构体数据
    CRTS_dma_iget(&slavePara, para, sizeof(MulAddPara), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);
    double beta = slavePara.beta_k;
    int cells = slavePara.cells;
    int id = CRTS_tid;

    int addr = slavePara.task[id].col_start;
    int len = slavePara.task[id].col_num;
    double *p = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *z = (double *)CRTS_pldm_malloc(len * sizeof(double));
    //接收数组数据
    CRTS_dma_iget(p, slavePara.p_k + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(z, slavePara.z_k1 + addr, len * sizeof(double), &DMARply);
    DMARplyCount += 2;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    //计算
    int i = 0;
    for (; i < len; i++) {
        p[i] = z[i] + beta * p[i];
    }
    //传回计算结果
    CRTS_dma_iput(slavePara.p_k + addr, p, len * sizeof(double), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    CRTS_pldm_free(p, len * sizeof(double));
    CRTS_pldm_free(z, len * sizeof(double));
}

// --------------------------------------------------------------------------------

void slave_Mul(MulPara *para) {
    MulPara slavePara;
    CRTS_dma_iget(&slavePara, para, sizeof(MulPara), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);
    int cells = slavePara.cells;
    int id = CRTS_tid;

    int addr = slavePara.task[id].col_start;
    int len = slavePara.task[id].col_num;
    double *M = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *r = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *z = (double *)CRTS_pldm_malloc(len * sizeof(double));

    CRTS_dma_iget(M, slavePara.m + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(r, slavePara.r_k1 + addr, len * sizeof(double), &DMARply);
    DMARplyCount += 2;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    for (int i = 0; i < len; ++i) {
        z[i] = M[i] * r[i];
    }

    CRTS_dma_iput(slavePara.z_k1 + addr, z, len * sizeof(double), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    CRTS_pldm_free(M, len * sizeof(double));
    CRTS_pldm_free(r, len * sizeof(double));
    CRTS_pldm_free(z, len * sizeof(double));
}

// --------------------------------------------------------------------------------

void slave_SubMul(SubMulPara *para) {
    SubMulPara slavePara;
    CRTS_dma_iget(&slavePara, para, sizeof(SubMulPara), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);
    int cells = slavePara.cells;
    int id = CRTS_tid;

    int addr = slavePara.task[id].col_start;
    int len = slavePara.task[id].col_num;
    double *M = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *r = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *g = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *z = (double *)CRTS_pldm_malloc(len * sizeof(double));

    CRTS_dma_iget(r, slavePara.r_k1 + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(g, slavePara.g + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(M, slavePara.m + addr, len * sizeof(double), &DMARply);
    DMARplyCount += 3;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    for (int i = 0; i < len; ++i) {
        z[i] = (r[i] - g[i]) * M[i];
    }

    CRTS_dma_iput(slavePara.z_k1 + addr, z, len * sizeof(double), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    CRTS_pldm_free(r, len * sizeof(double));
    CRTS_pldm_free(g, len * sizeof(double));
    CRTS_pldm_free(M, len * sizeof(double));
    CRTS_pldm_free(z, len * sizeof(double));
}

// --------------------------------------------------------------------------------

void slave_MulSub(MulSubPara *para) {
    MulSubPara slavePara;
    CRTS_dma_iget(&slavePara, para, sizeof(MulSubPara), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);
    int cells = slavePara.cells;
    int id = CRTS_tid;

    int addr = slavePara.task[id].col_start;
    int len = slavePara.task[id].col_num;
    double *g = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *M = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *z = (double *)CRTS_pldm_malloc(len * sizeof(double));

    CRTS_dma_iget(g, slavePara.g + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(z, slavePara.z_k1 + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(M, slavePara.m + addr, len * sizeof(double), &DMARply);
    DMARplyCount += 3;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    for (int i = 0; i < len; ++i) {
        g[i] = g[i] - z[i] * M[i];
    }

    CRTS_dma_iput(slavePara.g + addr, g, len * sizeof(double), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    CRTS_pldm_free(g, len * sizeof(double));
    CRTS_pldm_free(M, len * sizeof(double));
    CRTS_pldm_free(z, len * sizeof(double));
}

// --------------------------------------------------------------------------------

void slave_Reduce(ReducePara *para) {
    ReducePara slavePara;
    CRTS_dma_iget(&slavePara, para, sizeof(ReducePara), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);
    int cells = slavePara.cells;
    int id = CRTS_tid;

    int addr = slavePara.task[id].col_start;
    int len = slavePara.task[id].col_num;

    double *r = (double *)CRTS_pldm_malloc(len * sizeof(double));
    CRTS_dma_iget(r, slavePara.r_k1 + addr, len * sizeof(double), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    double local_sum = 0, sum = 0;
    for (int i = 0; i < len; ++i) {
        union {
            double f64;
            int64_t i64;
        } u;
        u.f64 = r[i];
        u.i64 &= 0x7fffffffffffffff;
        local_sum += u.f64;
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, &reducebuf, 64);
    if (CRTS_tid == 0) {
        int result = 1;
        CRTS_dma_put(slavePara.result, &result, sizeof(int));
        CRTS_dma_put(slavePara.residual, &sum, sizeof(double));
    }
    CRTS_pldm_free(r, len * sizeof(double));
}

// --------------------------------------------------------------------------------

void slave_MulReduceZR(MulReduceZRPara *para) {
    MulReduceZRPara slavePara;
    CRTS_dma_iget(&slavePara, para, sizeof(MulReduceZRPara), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);
    int cells = slavePara.cells;
    int id = CRTS_tid;

    int addr = slavePara.task[id].col_start;
    int len = slavePara.task[id].col_num;

    double *r = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *z = (double *)CRTS_pldm_malloc(len * sizeof(double));
    CRTS_dma_iget(r, slavePara.r_k1 + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(z, slavePara.z_k1 + addr, len * sizeof(double), &DMARply);
    DMARplyCount += 2;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    double local_sum = 0, sum = 0;
    for (int i = 0; i < len; ++i) {
        local_sum += r[i] * z[i];
    }

    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, &reducebuf, 64);
    if (CRTS_tid == 0)
        CRTS_dma_put(slavePara.result, &sum, sizeof(double));
    CRTS_pldm_free(r, len * sizeof(double));
    CRTS_pldm_free(z, len * sizeof(double));
}

// --------------------------------------------------------------------------------

void slave_MulReducepAx(MulReducepAxPara *para) {
    MulReducepAxPara slavePara;
    CRTS_dma_iget(&slavePara, para, sizeof(MulReducepAxPara), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);
    int cells = slavePara.cells;
    int id = CRTS_tid;

    int addr = slavePara.task[id].col_start;
    int len = slavePara.task[id].col_num;

    double *p = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *Ax = (double *)CRTS_pldm_malloc(len * sizeof(double));
    CRTS_dma_iget(p, slavePara.p_k + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(Ax, slavePara.Ax + addr, len * sizeof(double), &DMARply);
    DMARplyCount += 2;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    double local_sum = 0, sum = 0;
    for (int i = 0; i < len; ++i) {
        local_sum += p[i] * Ax[i];
    }

    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, &reducebuf, 64);
    if (CRTS_tid == 0)
        CRTS_dma_put(slavePara.result, &sum, sizeof(double));
    CRTS_pldm_free(p, len * sizeof(double));
    CRTS_pldm_free(Ax, len * sizeof(double));
}

// --------------------------------------------------------------------------------

void slave_Updatexr(UpdatexrPara *para) {
    UpdatexrPara slavePara;
    CRTS_dma_iget(&slavePara, para, sizeof(UpdatexrPara), &DMARply);
    DMARplyCount++;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);
    int cells = slavePara.cells;
    double alpha = slavePara.alpha;
    int id = CRTS_tid;

    int addr = slavePara.task[id].col_start;
    int len = slavePara.task[id].col_num;

    double *x = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *r = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *p = (double *)CRTS_pldm_malloc(len * sizeof(double));
    double *Ax = (double *)CRTS_pldm_malloc(len * sizeof(double));
    CRTS_dma_iget(x, slavePara.x + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(r, slavePara.r + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(p, slavePara.p + addr, len * sizeof(double), &DMARply);
    CRTS_dma_iget(Ax, slavePara.Ax + addr, len * sizeof(double), &DMARply);
    DMARplyCount += 4;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    for (int i = 0; i < len; ++i) {
        x[i] = x[i] + alpha * p[i];
        r[i] = r[i] - alpha * Ax[i];
    }

    CRTS_dma_iput(slavePara.x + addr, x, len * sizeof(double), &DMARply);
    CRTS_dma_iput(slavePara.r + addr, r, len * sizeof(double), &DMARply);
    DMARplyCount += 2;
    CRTS_dma_wait_value(&DMARply, DMARplyCount);

    CRTS_pldm_free(x, len * sizeof(double));
    CRTS_pldm_free(r, len * sizeof(double));
    CRTS_pldm_free(p, len * sizeof(double));
    CRTS_pldm_free(Ax, len * sizeof(double));
}
