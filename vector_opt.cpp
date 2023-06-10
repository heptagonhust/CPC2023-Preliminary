#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <crts.h>

#include "pcg_def.h"
#include "para_def.h"
#include "vector_opt.h"
#include "spmv_opt.h"

void pcg_init_precondition_csr_opt(const CsrMatrix &csr_matrix, Precondition &pre) {
    for(int i = 0 ; i < csr_matrix.rows; i++) {
        for(int j = csr_matrix.row_off[i]; j < csr_matrix.row_off[i+1]; j++){
            // get diagonal matrix
            if(csr_matrix.cols[j] == i) {
                pre.pre_mat_val[j] = 0.;	 
                pre.preD[i] = 1.0/csr_matrix.data[j];
            } else {
                pre.pre_mat_val[j] = csr_matrix.data[j];
            }
        }
    }
}

void pcg_precondition_csr_opt(const CsrMatrix &csr_matrix, const Precondition &pre, double *rAPtr, double *wAPtr) {
//! ---------------------------- START -----------------------------------------------
    MulPara para;
    para.m = pre.preD;
    para.r_k1 = rAPtr;
    para.z_k1 = wAPtr;
    para.cells = csr_matrix.rows;
    athread_spawn(Mul, &para);
    athread_join();
    // v_dot_product(csr_matrix.rows, pre.preD, rAPtr, wAPtr);
//! ------------------------------ END -----------------------------------------------
    double* gAPtr = (double*)malloc(csr_matrix.rows*sizeof(double));
    memset(gAPtr, 0, csr_matrix.rows*sizeof(double));
    for(int deg = 1; deg < 2; deg++) {
        // gAPtr = wAptr * pre.pre_mat_val; vec[rows] = matrix * vec[rows]
        csr_precondition_spmv_opt(csr_matrix, wAPtr, pre.pre_mat_val, gAPtr);
//! ---------------------------- START -----------------------------------------------
        SubMulPara para_sub;
        para_sub.r_k1 = rAPtr;
        para_sub.g = gAPtr;
        para_sub.m = pre.preD;
        para_sub.z_k1 = wAPtr;
        para_sub.cells = csr_matrix.rows;
        athread_spawn(SubMul, &para);
        athread_join();
        // v_sub_dot_product(csr_matrix.rows, rAPtr, gAPtr, pre.preD, wAPtr);
//! ------------------------------ END -----------------------------------------------
        memset(gAPtr, 0, csr_matrix.rows*sizeof(double));
    }
    free(gAPtr);
}

//! res += fasb(r[i])
double pcg_gsumMag_opt(double *r, int size, double normfactor, double tolerance) {
    int result;
    double residual;
    ReducePara para;
    para.cells = size;
    para.normfactor = normfactor;
    para.tolerance = tolerance;
    para.r_k1 = r;
    para.result = &result;
    para.residual = &residual;
    athread_spawn(Reduce, &para);
    athread_join();
    return residual;
}


// multiply and reduce : vector inner product
// 逐元素与规约操作，需要核间通信
// 存在不同的参数调用，分离pcg_gsumProd成为两个函数分别进行调用
//! result += z[i] * r[i]
double pcg_gsumProd_opt_zr(double *z, double *r, int size) {
    MulReduceZRPara para;
    double sum = .0;
    para.cells = size;
    para.z_k1 = z;
    para.r_k1 = r;
    para.result = &sum;
    athread_spawn(MulReduceZR, &para);
    athread_join();
    return sum;
}

//! result += p[i] * Ax[i]
double pcg_gsumProd_opt_pAx(double *p, double *Ax, int size) {
    MulReducepAxPara para;
    double sum = .0;
    para.cells = size;
    para.p_k = p;
    para.Ax = Ax;
    para.result = &sum;
    athread_spawn(MulReducepAx, &para);
    athread_join();
    return sum;
}

//! x = x + alpha * p
//! r = r - alpha * Ax
void pcg_update_xr_opt(double *x, double *r, double *p, double *Ax, double alpha, int cells) {
    UpdatexrPara para;
    para.cells = cells;
    para.alpha = alpha;
    para.x = x;
    para.r = r;
    para.p = p;
    para.Ax = Ax;
    athread_spawn(Updatexr, &para);
    athread_join();
}