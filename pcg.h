#ifndef _PCG_H_
#define _PCG_H_

#include "pcg_def.h"
#include "vector.h"

void read_mesh();

void init_p_equation(
    LduMatrix &ldu_matrix,
    double *&source,
    double *&psi,
    int mesh);
void free_p_equation(LduMatrix &ldu_matrix, double *source, double *psi);
int check_result(const PCGReturn &pcg_return, double run_time, int mesh);
void write_result(const PCGReturn &pcg_return, int mesh);

// PCG
PCGReturn pcg_solve(
    const LduMatrix &ldu_matrix,
    double *source,
    double *psi,
    int maxIter,
    double tolerance,
    double normfactor);
void pcg_init_precondition_csc(
    const CsrMatrix &csr_matrix,
    Precondition &pre,
    double *M);
void free_pcg(PCG &pcg);
void free_csr_matrix(CsrMatrix &csr_matrix);
void free_precondition(Precondition &pre);
void pcg_init_precondition_csc(
    const CscMatrix &csr_matrix,
    Precondition &pre,
    double *M);
void pcg_precondition_csc_opt(
    SpmvPara &para_Az,
    const Precondition &pre,
    double *rAPtr,
    double *wAPtr,
    double *M,
    Slave_task *ntask,
    int cells);
#endif
