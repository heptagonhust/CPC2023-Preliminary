#pragma once
#include "pcg_def.h"

void csr_spmv_opt(const CsrMatrix &csr_matrix, double *x, double *b);
void csr_precondition_spmv_opt(
    const CsrMatrix &csr_matrix,
    double *vec,
    double *val,
    double *result);