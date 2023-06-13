#include "vector_master.h"

#include <crts.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cstring>

#include "pcg.h"
#include "pcg_def.h"


//! res += fasb(r[i])
double pcg_gsumMag_opt(
    double *r,
    int size,
    double normfactor,
    double tolerance,
    Slave_task *ntask) {
    int result;
    double residual;
    ReducePara para;
    para.cells = size;
    para.normfactor = normfactor;
    para.tolerance = tolerance;
    para.r_k1 = r;
    para.result = &result;
    para.residual = &residual;
    memcpy(&para.task, ntask, 64 * sizeof(Slave_task));
    athread_spawn(slave_Reduce, &para);
    athread_join();
    // printf("vector_opt, 77: residual: %f\n", residual);
    return residual;
}

// multiply and reduce : vector inner product
// 逐元素与规约操作，需要核间通信
// 存在不同的参数调用，分离pcg_gsumProd成为两个函数分别进行调用
//! result += z[i] * r[i]
double pcg_gsumProd_opt_zr(double *z, double *r, int size, Slave_task *ntask) {
    MulReduceZRPara para;
    double sum = .0;
    para.cells = size;
    para.z_k1 = z;
    para.r_k1 = r;
    para.result = &sum;
    memcpy(&para.task, ntask, 64 * sizeof(Slave_task));
    athread_spawn(slave_MulReduceZR, &para);
    athread_join();
    return sum;
}

//! result += p[i] * Ax[i]
double
pcg_gsumProd_opt_pAx(double *p, double *Ax, int size, Slave_task *ntask) {
    MulReducepAxPara para;
    double sum = .0;
    para.cells = size;
    para.p_k = p;
    para.Ax = Ax;
    para.result = &sum;
    memcpy(&para.task, ntask, 64 * sizeof(Slave_task));
    athread_spawn(slave_MulReducepAx, &para);
    athread_join();
    return sum;
}

void v_dot_product_opt(
    const int nCells,
    double *m,
    double *r_k1,
    double *z_k1,
    Slave_task *ntask) {
    MulPara para;
    para.cells = nCells;
    para.m = m;
    para.r_k1 = r_k1;
    para.z_k1 = z_k1;
    memcpy(&para.task, ntask, 64 * sizeof(Slave_task));
    athread_spawn(slave_Mul, &para);
    athread_join();
}

void v_sub_dot_product_opt(
    const int nCells,
    double *r_k1,
    double *g,
    double *m,
    double *z_k1,
    Slave_task *ntask) {
    SubMulPara para;
    para.cells = nCells;
    para.r_k1 = r_k1;
    para.g = g;
    para.m = m;
    para.z_k1 = z_k1;
    memcpy(&para.task, ntask, 64 * sizeof(Slave_task));
    athread_spawn(slave_SubMul, &para);
    athread_join();
}

//! x = x + alpha * p
//! r = r - alpha * Ax
void pcg_update_xr_opt(
    double *x,
    double *r,
    double *p,
    double *Ax,
    double alpha,
    int cells,
    Slave_task *ntask) {
    UpdatexrPara para;
    para.cells = cells;
    para.alpha = alpha;
    para.x = x;
    para.r = r;
    para.p = p;
    para.Ax = Ax;
    memcpy(&para.task, ntask, 64 * sizeof(Slave_task));
    athread_spawn(slave_Updatexr, &para);
    athread_join();
}

void pcg_update_p_opt(
    double *p,
    double *z,
    double beta,
    int cells,
    Slave_task *ntask) {
    MulAddPara para;
    para.p_k = p;
    para.z_k1 = z;
    para.beta_k = beta;
    para.cells = cells, memcpy(&para.task, ntask, 64 * sizeof(Slave_task));
    athread_spawn(slave_MulAdd, &para);
    athread_join();
}
