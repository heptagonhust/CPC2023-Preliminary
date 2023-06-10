#pragma once

#include "pcg_def.h"
#include "para_def.h"

extern "C" void MulAdd(MulAddPara *para);
extern "C" void Mul(MulPara *para);
extern "C" void SubMul(SubMulPara *para);
extern "C" void Reduce(ReducePara *para);
extern "C" void MulReduceZR(MulReduceZRPara *para);
extern "C" void MulReducepAx(MulReducepAxPara *para);
extern "C" void Updatexr(UpdatexrPara *para);

void pcg_init_precondition_csr_opt(
    const CsrMatrix &csr_matrix,
    Precondition &pre);
void pcg_precondition_csr_opt(
    const CsrMatrix &csr_matrix,
    const Precondition &pre,
    double *rAPtr,
    double *wAPtr);
double pcg_gsumMag_opt(double *r, int size, double normfactor, double tolerance);
double pcg_gsumProd_opt_zr(double *z, double *r, int size);
double pcg_gsumProd_opt_pAx(double *z, double *r, int size);
void pcg_update_xr_opt(
    double *x,
    double *r,
    double *p,
    double *Ax,
    double alpha,
    int cells);

void csr_spmv_opt(const CsrMatrix &csr_matrix, double *x, double *b);
void csr_precondition_spmv_opt(
    const CsrMatrix &csr_matrix,
    double *vec,
    double *val,
    double *result);