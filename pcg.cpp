#include "pcg.h"

#include <crts.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coo_matrix.h"

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

    // Precondition pre;
    // pre.preD = (double *)malloc(cells * sizeof(double));
    // pre.pre_mat_val = (double *)malloc((cells + faces * 2) * sizeof(double));
    double *preD = (double *)malloc(cells * sizeof(double));

    double *M = (double *)malloc(cells * sizeof(double));

    CRTS_init();

    Slave_task ntask[64];
    int len = cells / 64;
    int rest = cells % 64;
    for (int i = 0; i < 64; ++i) {
        if (i < rest) {
            ntask[i].col_num = len + 1;
            ntask[i].col_start = i * (len + 1);
        } else {
            ntask[i].col_num = len;
            ntask[i].col_start = i * len + rest;
        }
    }

    // data structure transform
    struct timeval tv;
	gettimeofday(&tv,NULL);
    double start = (double)(tv.tv_sec)+(double)(tv.tv_usec)*1e-6;
    SplitedCooMatrix splited_matrix = ldu_to_splited_coo(ldu_matrix, SLAVE_CORE_NUM);
    gettimeofday(&tv,NULL);
    double end = (double)(tv.tv_sec)+(double)(tv.tv_usec)*1e-6;
    INFO("split time: %.4lfs\n", end - start);

    // Spmv parameter generation
    SpmvPara para_Ax, para_Ap, para_Az;
    spmv_para_from_splited_coo_matrix(&splited_matrix, &para_Ap, pcg.p, cells, cells);
    spmv_para_from_splited_coo_matrix(&splited_matrix, &para_Az, pcg.z, cells, cells);
    spmv_para_from_splited_coo_matrix(&splited_matrix, &para_Ax, x, cells, cells);

    pcg_init_precondition_ldu(ldu_matrix, preD, M);

    // AX = A * X
    coo_spmv(&para_Ax, pcg.Ax);
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
                pcg_precondition_coo_opt(para_Az, preD, pcg.r, pcg.z, M, ntask, cells);
                // tol_0= swap(r) * z
                pcg.sumprod = pcg_gsumProd_opt_zr(pcg.r, pcg.z, cells, ntask);
                // iter ==0 ; p = z
                memcpy(pcg.p, pcg.z, cells * sizeof(double));
            } else {
                pcg.sumprod_old = pcg.sumprod;
                // z = M(-1) * r
                pcg_precondition_coo_opt(para_Az, preD, pcg.r, pcg.z, M, ntask, cells);
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
            coo_spmv(&para_Ap, pcg.Ax);

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

    free(preD);
    free_pcg(pcg);

    spmv_para_free(&para_Ax);
    spmv_para_free(&para_Ap);
    spmv_para_free(&para_Az);

    PCGReturn pcg_return;
    pcg_return.residual = pcg.residual;
    pcg_return.iter = iter;
    return pcg_return;
}


// // diagonal precondition, get matrix M^(-1) (diagonal matrix)
// // pre_mat_val: 非对角元     : csr_matrix中元素
// //              对角元素     : 0
// // preD       : csr_matrix中对角元素的倒数
// void pcg_init_precondition_csc(const CscMatrix &csc_matrix, Precondition &pre, double *M) {
//     for (int i = 0; i < csc_matrix.cols; i++) {
//         for (int j = csc_matrix.col_off[i]; j < csc_matrix.col_off[i + 1];
//              j++) {
//             if (csc_matrix.rows[j] == i) {
//                 pre.preD[i] = 1.0 / csc_matrix.data[j];
//                 M[i] = csc_matrix.data[j];
//             }
//         }
//     }
// }

void pcg_init_precondition_ldu(
    const LduMatrix &ldu_matrix,
    double *preD,
    double *M) {
    for (int i = 0; i < ldu_matrix.cells; i++) {
        preD[i] = 1.0 / ldu_matrix.diag[i];
        M[i] = ldu_matrix.diag[i];
    }
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
    memset(gAPtr, 0, cells * sizeof(double));
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
    free(pcg.z);
    free(pcg.p);
    free(pcg.Ax);
}

void free_precondition(Precondition &pre) {
    free(pre.preD);
    free(pre.pre_mat_val);
}
