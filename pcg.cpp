#include "pcg.h"

#include <crts.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "csc_matrix.cpp"
#include "csr_matrix.cpp"
#include "vector_utils.cpp"

// 示例
typedef struct {
    double *p;
    double *z;
    double beta;
    int cells;
} Para;
extern "C" void slave_example(Para *para);

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

    // format transform
    CsrMatrix csr_matrix;
    ldu_to_csr(ldu_matrix, csr_matrix);

    pcg_init_precondition_csr(csr_matrix, pre);

    // AX = A * X
    csr_spmv(csr_matrix, x, pcg.Ax);
    // r = b - A * x
    for (int i = 0; i < cells; i++) {
        pcg.r[i] = source[i] - pcg.Ax[i];
    }
    // calculate residual, scale
    pcg.residual = pcg_gsumMag(pcg.r, cells);
    double init_residual = pcg.residual;

    if (fabs(pcg.residual / normfactor) > tolerance) {
        do {
            if (iter == 0) {
                // z = M(-1) * r
                // M: diagonal matrix of csr matrix A : diagonal preprocess
                pcg_precondition_csr(csr_matrix, pre, pcg.r, pcg.z);
                // tol_0= swap(r) * z
                pcg.sumprod = pcg_gsumProd(pcg.r, pcg.z, cells);
                // iter ==0 ; p = z
                memcpy(pcg.p, pcg.z, cells * sizeof(double));
            } else {
                pcg.sumprod_old = pcg.sumprod;
                // z = M(-1) * r
                pcg_precondition_csr(csr_matrix, pre, pcg.r, pcg.z);
                // tol_0= swap(r) * z
                pcg.sumprod = pcg_gsumProd(pcg.r, pcg.z, cells);
                // beta = tol_1 / tol_0
                // p = z + beta * p
                pcg.beta = pcg.sumprod / pcg.sumprod_old;

                // 未优化代码段
                /*for(int i = 0; i < cells; i++) {
                    pcg.p[i] = pcg.z[i] + pcg.beta * pcg.p[i];
                }*/

                // == 优化示例代码段 ==
                static int isInit = 0;
                if (isInit == 0) {
                    // 从核初始化
                    CRTS_init();
                    isInit = 1;
                }
                // 参数定义并赋值
                Para para;
                para.p = pcg.p;
                para.z = pcg.z;
                para.beta = pcg.beta;
                para.cells = cells;
                // 启动从核
                athread_spawn(slave_example, &para);
                // 等待从核线程组终止
                athread_join();
                // == 优化示例代码段 ==
            }

            // Ax = A * p
            csr_spmv(csr_matrix, pcg.p, pcg.Ax);

            // alpha = tol_0 / tol_1 = (swap(r) * z) / ( swap(p) * A * p)
            pcg.alpha = pcg.sumprod / pcg_gsumProd(pcg.p, pcg.Ax, cells);

            // x = x + alpha * p
            // r = r - alpha * Ax
            for (int i = 0; i < cells; i++) {
                x[i] = x[i] + pcg.alpha * pcg.p[i];
                pcg.r[i] = pcg.r[i] - pcg.alpha * pcg.Ax[i];
            }

            // tol_1 = swap(z) * r
            pcg.residual = pcg_gsumMag(pcg.r, cells);
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


// reduce
// 规约操作，需要核间通信
double pcg_gsumMag(double *r, int size) {
    double ret = .0;
    for (int i = 0; i < size; i++) {
        ret += fabs(r[i]);
    }
    return ret;
}

// multiply and reduce : vector inner product
// 逐元素与规约操作，需要核间通信
double pcg_gsumProd(double *z, double *r, int size) {
    double ret = .0;
    for (int i = 0; i < size; i++) {
        ret += z[i] * r[i];
    }
    return ret;
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
