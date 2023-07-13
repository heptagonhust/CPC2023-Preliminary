#include "pcg_opt.h"

#include <crts.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coo_matrix.h"

#include "vector.h"
#include "spmv.h"
#include "slave_def.h"

#define SLAVE_CORE_NUM 64

extern "C" void slave_MainLoop(MainLoopPara *para_mem);

// ldu_matrix: matrix A
// source: vector b
PCGReturn pcg_solve(
    const LduMatrix &ldu_matrix,
    double *source,
    double *x,
    int maxIter,
    double tolerance,
    double normfactor) {
    int iter = 0;
    // cells: matrix rows
    int cells = ldu_matrix.cells;
    int faces = ldu_matrix.faces;
    double init_residual;

    PCG pcg;
    pcg.r = (double *)malloc(cells * sizeof(double));
    pcg.p = (double *)malloc(cells * sizeof(double));
    pcg.x = x;
    pcg.source = source;

    double *preD = (double *)malloc(cells * sizeof(double));
    double *M = ldu_matrix.diag;

    CRTS_init();

    // data structure transform
    struct timeval tv;
	gettimeofday(&tv,NULL);
    double start = (double)(tv.tv_sec)+(double)(tv.tv_usec)*1e-6;
    SplitedCooMatrix splited_matrix = ldu_to_splited_coo(ldu_matrix, SLAVE_CORE_NUM);
    gettimeofday(&tv,NULL);
    double end = (double)(tv.tv_sec)+(double)(tv.tv_usec)*1e-6;
    // INFO("split time: %.4lfs\n", end - start);

    Slave_task ntask[SLAVE_CORE_NUM];
    for (int i = 0; i < SLAVE_CORE_NUM; ++i) {
        ntask[i].col_num = splited_matrix.chunk_ranges[i].row_num;
        ntask[i].col_start = splited_matrix.chunk_ranges[i].row_begin;
    }

    // Spmv parameter generation
    SpmvPara para_A;
    spmv_para_from_splited_coo_matrix(&splited_matrix, &para_A, cells, cells);

    memcpy(pcg.p, pcg.x, cells * sizeof(double));
    memcpy(pcg.r, source, cells * sizeof(double));

    MainLoopPara para;
    para.maxIter = maxIter;
    para.normfactor = normfactor;
    para.tolerance = tolerance;
    para.p = pcg.p;
    para.r = pcg.r;
    para.x = pcg.x;
    para.M = M;
    para.spmv_para = &para_A;
    memcpy(&para.ntask, &ntask, SLAVE_CORE_NUM * sizeof(Slave_task));
    para.init_residual = &init_residual;
    para.final_residual = &pcg.residual;
    para.iter = &iter;
    athread_spawn(slave_MainLoop, &para);
    athread_join();

    INFO(
        "PCG: init residual = %e, final residual = %e, iterations: %d\n",
        init_residual,
        pcg.residual,
        iter);

    free(preD);
    free_pcg(pcg);

    PCGReturn pcg_return;
    pcg_return.residual = pcg.residual;
    pcg_return.iter = iter;
    return pcg_return;
}

void pcg_precondition_coo_opt(
    SpmvPara &para_Az,
    double *preD,
    double *rAPtr,
    double *wAPtr,
    double *M,
    Slave_task *ntask,
    int cells) {
    v_dot_product_opt(cells, preD, rAPtr, wAPtr, ntask);
    double *gAPtr = (double *)malloc(cells * sizeof(double));
    for (int deg = 1; deg < 2; deg++) {
        coo_spmv(&para_Az, gAPtr);
        precond_update_g_opt(gAPtr, wAPtr, M, cells, ntask);
        v_sub_dot_product_opt(
            cells,
            rAPtr,
            gAPtr,
            preD,
            wAPtr,
            ntask);
        memset(gAPtr, 0, cells * sizeof(double));
    }
    free(gAPtr);
}


void free_pcg(PCG &pcg) {
    free(pcg.r);
    free(pcg.p);
}

void free_precondition(Precondition &pre) {
    free(pre.preD);
    free(pre.pre_mat_val);
}
