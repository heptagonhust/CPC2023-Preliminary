#pragma once

#include "pcg_def.h"
#include "slave_def.h"

void pcg_init_precondition_csr_opt(
    const CsrMatrix &csr_matrix,
    Precondition &pre,
    Slave_task *ntask);
void pcg_precondition_csr_opt(
    const CsrMatrix &csr_matrix,
    const Precondition &pre,
    double *rAPtr,
    double *wAPtr,
    Slave_task *ntask);
double pcg_gsumMag_opt(
    double *r,
    int size,
    double normfactor,
    double tolerance,
    Slave_task *ntask);
double pcg_gsumProd_opt_zr(double *z, double *r, int size, Slave_task *ntask);
double pcg_gsumProd_opt_pAx(double *p, double *Ax, int size, Slave_task *ntask);
void v_sub_dot_product_opt(
    int nCells,
    double *r_k1,
    double *g,
    double *m,
    double *z_k1,
    Slave_task *ntask);
void v_dot_product_opt(
    int nCells,
    double *m,
    double *r_k1,
    double *z_k1,
    Slave_task *ntask);
void pcg_update_xr_opt(
    double *x,
    double *r,
    double *p,
    double *Ax,
    double alpha,
    int cells,
    Slave_task *ntask);
void pcg_update_p_opt(
    double *p,
    double *z,
    double beta,
    int cells,
    Slave_task *ntask);
