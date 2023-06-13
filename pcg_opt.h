#pragma once

#include "pcg_def.h"
#include "vector/vector.h"


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