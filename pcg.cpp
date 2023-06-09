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
typedef struct
{
    double *p;
    double *z;
    double beta;
    int cells;
} Para;

typedef struct
{
    // TODO
} CscChunkPara;

extern "C" void slave_example(Para *para);
extern "C" void csc_spmv_slave(CscChunkPara *chunk);

// ldu_matrix: matrix A
// source: vector b
PCGReturn pcg_solve(
    const LduMatrix &ldu_matrix,
    double *source,
    double *x,
    int maxIter,
    double tolerance,
    double normfactor)
{
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
    for (int i = 0; i < cells; i++)
    {
        pcg.r[i] = source[i] - pcg.Ax[i];
    }
    // calculate residual, scale
    pcg.residual = pcg_gsumMag(pcg.r, cells);
    double init_residual = pcg.residual;

    if (fabs(pcg.residual / normfactor) > tolerance)
    {
        do
        {
            if (iter == 0)
            {
                // z = M(-1) * r
                // M: diagonal matrix of csr matrix A : diagonal preprocess
                pcg_precondition_csr(csr_matrix, pre, pcg.r, pcg.z);
                // tol_0= swap(r) * z
                pcg.sumprod = pcg_gsumProd(pcg.r, pcg.z, cells);
                // iter ==0 ; p = z
                memcpy(pcg.p, pcg.z, cells * sizeof(double));
            }
            else
            {
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
                if (isInit == 0)
                {
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
            for (int i = 0; i < cells; i++)
            {
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

void ldu_to_csr(const LduMatrix &ldu_matrix, CsrMatrix &csr_matrix)
{
    csr_matrix.rows = ldu_matrix.cells;
    csr_matrix.data_size = 2 * ldu_matrix.faces + ldu_matrix.cells;
    csr_matrix.row_off = (int *)malloc((csr_matrix.rows + 1) * sizeof(int));
    csr_matrix.cols = (int *)malloc(csr_matrix.data_size * sizeof(int));
    csr_matrix.data = (double *)malloc(csr_matrix.data_size * sizeof(double));

    int row, col, offset;
    int *tmp = (int *)malloc((csr_matrix.rows + 1) * sizeof(int));

    csr_matrix.row_off[0] = 0;
    for (int i = 1; i < csr_matrix.rows + 1; i++)
        csr_matrix.row_off[i] = 1;

    for (int i = 0; i < ldu_matrix.faces; i++)
    {
        row = ldu_matrix.uPtr[i];
        col = ldu_matrix.lPtr[i];
        csr_matrix.row_off[row + 1]++;
        csr_matrix.row_off[col + 1]++;
    }

    for (int i = 0; i < ldu_matrix.cells; i++)
    {
        csr_matrix.row_off[i + 1] += csr_matrix.row_off[i];
    }

    memcpy(&tmp[0], &csr_matrix.row_off[0], (ldu_matrix.cells + 1) * sizeof(int));
    // lower
    for (int i = 0; i < ldu_matrix.faces; i++)
    {
        row = ldu_matrix.uPtr[i];
        col = ldu_matrix.lPtr[i];
        offset = tmp[row]++;
        csr_matrix.cols[offset] = col;
        csr_matrix.data[offset] = ldu_matrix.lower[i];
    }

    // diag
    for (int i = 0; i < ldu_matrix.cells; i++)
    {
        offset = tmp[i]++;
        csr_matrix.cols[offset] = i;
        csr_matrix.data[offset] = ldu_matrix.diag[i];
    }

    // upper
    for (int i = 0; i < ldu_matrix.faces; i++)
    {
        row = ldu_matrix.lPtr[i];
        col = ldu_matrix.uPtr[i];
        offset = tmp[row]++;
        csr_matrix.cols[offset] = col;
        csr_matrix.data[offset] = ldu_matrix.upper[i];
    }

    free(tmp);
}

// basic spmv, 需要负载均衡
void csr_spmv(const CsrMatrix &csr_matrix, double *vec, double *result)
{
    for (int i = 0; i < csr_matrix.rows; i++)
    {
        int start = csr_matrix.row_off[i];
        int num = csr_matrix.row_off[i + 1] - csr_matrix.row_off[i];
        double temp = 0;
        for (int j = 0; j < num; j++)
        {
            temp += vec[csr_matrix.cols[start + j]] * csr_matrix.data[start + j];
        }
        result[i] = temp;
    }
}

#define N_SLAVE_CORES 64

void csc_spmv(const CscMatrix &csc_matrix, double *vec, double *result)
{
    int chunk_size = csc_matrix.cols / N_SLAVE_CORES;
    // #ifndef NDEBUG
    //     printf("current slave cores: %d", CRTS_athread_get_max_threads());
    // #endif
    for (int i = 0; i < N_SLAVE_CORES; ++i)
    {
        CscChunkPara chunk_para;
        // TODO: 将 csc 分成 64 个 chunk
        CRTS_athread_create(i, csc_spmv_slave, &chunk_para);
    }

    // 主核与核组同步
    CRTS_sync_master_array();
}

void csr_precondition_spmv(const CsrMatrix &csr_matrix, double *vec, double *val, double *result)
{
    for (int i = 0; i < csr_matrix.rows; i++)
    {
        int start = csr_matrix.row_off[i];
        int num = csr_matrix.row_off[i + 1] - csr_matrix.row_off[i];
        double temp = 0;
        for (int j = 0; j < num; j++)
        {
            temp += vec[csr_matrix.cols[start + j]] * val[start + j];
        }
        result[i] = temp;
    }
}

void v_dot_product(const int nCells, const double *vec1, const double *vec2, double *result)
{
    for (int cell = 0; cell < nCells; cell++)
    {
        result[cell] = vec1[cell] * vec2[cell];
    }
}

void v_sub_dot_product(const int nCells, const double *sub, const double *subed, const double *vec, double *result)
{
    for (int cell = 0; cell < nCells; cell++)
    {
        result[cell] = (sub[cell] - subed[cell]) * vec[cell];
    }
}

// diagonal precondition, get matrix M^(-1) (diagonal matrix)
// pre_mat_val: 非对角元     : csr_matrix中元素
//              对角元素     : 0
// preD       : csr_matrix中对角元素的倒数
void pcg_init_precondition_csr(const CsrMatrix &csr_matrix, Precondition &pre)
{
    for (int i = 0; i < csr_matrix.rows; i++)
    {
        for (int j = csr_matrix.row_off[i]; j < csr_matrix.row_off[i + 1]; j++)
        {
            // get diagonal matrix
            if (csr_matrix.cols[j] == i)
            {
                pre.pre_mat_val[j] = 0.;
                pre.preD[i] = 1.0 / csr_matrix.data[j];
            }
            else
            {
                pre.pre_mat_val[j] = csr_matrix.data[j];
            }
        }
    }
}

// ? 存疑，循环用处?
void pcg_precondition_csr(const CsrMatrix &csr_matrix, const Precondition &pre, double *rAPtr, double *wAPtr)
{
    double *gAPtr = (double *)malloc(csr_matrix.rows * sizeof(double));
    v_dot_product(csr_matrix.rows, pre.preD, rAPtr, wAPtr);
    memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
    for (int deg = 1; deg < 2; deg++)
    {
        // gAPtr = wAptr * pre.pre_mat_val; vec[rows] = matrix * vec[rows]
        csr_precondition_spmv(csr_matrix, wAPtr, pre.pre_mat_val, gAPtr);
        v_sub_dot_product(csr_matrix.rows, rAPtr, gAPtr, pre.preD, wAPtr);
        memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
    }
    free(gAPtr);
}

// reduce
// 规约操作，需要核间通信
double pcg_gsumMag(double *r, int size)
{
    double ret = .0;
    for (int i = 0; i < size; i++)
    {
        ret += fabs(r[i]);
    }
    return ret;
}

// multiply and reduce : vector inner product
// 逐元素与规约操作，需要核间通信
double pcg_gsumProd(double *z, double *r, int size)
{
    double ret = .0;
    for (int i = 0; i < size; i++)
    {
        ret += z[i] * r[i];
    }
    return ret;
}

void free_pcg(PCG &pcg)
{
    free(pcg.r);
    free(pcg.z);
    free(pcg.p);
    free(pcg.Ax);
}

void free_precondition(Precondition &pre)
{
    free(pre.preD);
    free(pre.pre_mat_val);
}
