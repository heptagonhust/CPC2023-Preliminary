#include "pcg.h"

#include <crts.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "csc_matrix.cpp"
#include "vector_utils.cpp"

#include "vector.h"
#include "spmv.h"

#define SLAVE_CORE_NUM 64

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

    PCG pcg;
    pcg.r = (double *)malloc(cells * sizeof(double));
    pcg.z = (double *)malloc(cells * sizeof(double));
    pcg.p = (double *)malloc(cells * sizeof(double));
    pcg.Ax = (double *)malloc(cells * sizeof(double));
    pcg.x = x;
    pcg.source = source;

    Precondition pre;
    pre.preD = (double *)malloc(cells * sizeof(double));
    pre.pre_mat_val = (double *)malloc((cells + faces * 2) * sizeof(double));
    
    double *M = (double *)malloc(cells * sizeof(double));

    CRTS_init();

    // Slave_task ntask[64];
    // int len = cells / 64;
    // int rest = cells % 64;
    // for (int i = 0; i < 64; ++i) {
    //     if (i < rest) {
    //         ntask[i].col_num = len + 1;
    //         ntask[i].col_start = i * (len + 1);
    //     } else {
    //         ntask[i].col_num = len;
    //         ntask[i].col_start = i * len + rest;
    //     }
    // }

    // data structure transform
    CscMatrix csc_matrix;
    ldu_to_csc(ldu_matrix, csc_matrix);
    SplitedCscMatrix splited_matrix = split_csc_matrix(csc_matrix, SLAVE_CORE_NUM, NULL);
    // Spmv parameter generation
    SpmvPara para_Ax, para_Ap, para_Az;
    spmv_para_from_splited_csc_matrix(splited_matrix, para_Ap, pcg.p, cells, cells);
    spmv_para_from_splited_csc_matrix(splited_matrix, para_Az, pcg.z, cells, cells);
    spmv_para_from_splited_csc_matrix(splited_matrix, para_Ax, x, cells, cells);

    pcg_init_precondition_csc(csc_matrix, pre, M);

    // AX = A * X
    csc_spmv(para_Ax, pcg.Ax);
    // r = b - A * x
    for (int i = 0; i < cells; i++) {
        pcg.r[i] = source[i] - pcg.Ax[i];
    }
    // calculate residual, scale
    pcg.residual = pcg_gsumMag_opt(pcg.r, cells, normfactor, tolerance, ntask);

    double init_residual = pcg.residual;

    if (fabs(pcg.residual / normfactor) > tolerance) {
        do {
            if (iter == 0) {
                // z = M(-1) * r
                // M: diagonal matrix of csr matrix A : diagonal preprocess
                pcg_precondition_csr_opt(csr_matrix, pre, pcg.r, pcg.z, M, ntask);
                // tol_0= swap(r) * z
                pcg.sumprod = pcg_gsumProd_opt_zr(pcg.r, pcg.z, cells, ntask);
                // iter ==0 ; p = z
                memcpy(pcg.p, pcg.z, cells * sizeof(double));
            } else {
                pcg.sumprod_old = pcg.sumprod;
                // z = M(-1) * r
                pcg_precondition_csr_opt(csr_matrix, pre, pcg.r, pcg.z, M, ntask);
                // tol_0= swap(r) * z
                pcg.sumprod = pcg_gsumProd_opt_zr(pcg.r, pcg.z, cells, ntask);
                // beta = tol_1 / tol_0
                // p = z + beta * p
                pcg.beta = pcg.sumprod / pcg.sumprod_old;

                // 未优化代码段
                /*for(int i = 0; i < cells; i++) {
                    pcg.p[i] = pcg.z[i] + pcg.beta * pcg.p[i];
                }*/


                pcg_update_p_opt(pcg.p, pcg.z, pcg.beta, cells, ntask);
            }

            // Ax = A * p
            csc_spmv(para_Ap, pcg.Ax);

            // alpha = tol_0 / tol_1 = (swap(r) * z) / ( swap(p) * A * p)
            pcg.alpha =
                pcg.sumprod / pcg_gsumProd_opt_pAx(pcg.p, pcg.Ax, cells, ntask);

            // x = x + alpha * p
            // r = r - alpha * Ax
            pcg_update_xr_opt(x, pcg.r, pcg.p, pcg.Ax, pcg.alpha, cells, ntask);

            // tol_1 = swap(z) * r
            pcg.residual =
                pcg_gsumMag_opt(pcg.r, cells, normfactor, tolerance, ntask);
        } while (++iter < maxIter && (pcg.residual / normfactor) >= tolerance);
    }

    INFO(
        "PCG: init residual = %e, final residual = %e, iterations: %d\n",
        init_residual,
        pcg.residual,
        iter);

    free_pcg(pcg);
    free_csr_matrix(csr_matrix);
    free_precondition(pre);

    PCGReturn pcg_return;
    pcg_return.residual = pcg.residual;
    pcg_return.iter = iter;
    return pcg_return;
}


// diagonal precondition, get matrix M^(-1) (diagonal matrix)
// pre_mat_val: 非对角元     : csr_matrix中元素
//              对角元素     : 0
// preD       : csr_matrix中对角元素的倒数
void pcg_init_precondition_csc(const CscMatrix &csc_matrix, Precondition &pre, double *M) {
    for (int i = 0; i < csc_matrix.cols; i++) {
        for (int j = csr_matrix.col_off[i]; j < csc_matrix.col_off[i + 1];
             j++) {
            if (csc_matrix.rows[j] == i) {
                pre.preD[i] = 1.0 / csr_matrix.data[j];
                M[i] = csr_matrix.data[j];
            }
        }
    }
}

void pcg_precondition_csr_opt(
    const SpmvPara &para_Az,
    const Precondition &pre,
    double *rAPtr,
    double *wAPtr,
    double *M,
    Slave_task *ntask) {
    v_dot_product_opt(csr_matrix.rows, pre.preD, rAPtr, wAPtr, ntask);
    double *gAPtr = (double *)malloc(csr_matrix.rows * sizeof(double));
    memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
    for (int deg = 1; deg < 2; deg++) {
        csc_spmv(para_Az, gAPtr);
        precond_update_g_opt(gAPtr, wAPtr, M, csr_matrix.rows, ntask);
        v_sub_dot_product_opt(
            csr_matrix.rows,
            rAPtr,
            gAPtr,
            pre.preD,
            wAPtr,
            ntask);
        memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
    }
    free(gAPtr);
}


void free_pcg(PCG &pcg) {
    free(pcg.r);
    free(pcg.z);
    free(pcg.p);
    free(pcg.Ax);
}

void free_precondition(Precondition &pre) {
    free(pre.preD);
    free(pre.pre_mat_val);
}
