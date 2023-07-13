#pragma once

#include "pcg.h"
#include "vector.h"
#include "pcg_def_opt.h"

void pcg_init_precondition_ldu(
    const LduMatrix &ldu_matrix,
    double *preD,
    double *M);
void pcg_precondition_coo_opt(
    SpmvPara &para_Az,
    double *preD,
    double *rAPtr,
    double *wAPtr,
    double *M,
    Slave_task *ntask,
    int cells);