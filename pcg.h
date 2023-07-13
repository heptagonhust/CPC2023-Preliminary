#ifndef _PCG_H_
#define _PCG_H_

#include "pcg_def.h"

void read_mesh();

void init_p_equation(LduMatrix &ldu_matrix, double* &source, double* &psi, int mesh);
void free_p_equation(LduMatrix &ldu_matrix, double *source, double *psi);
int check_result(const PCGReturn &pcg_return, double run_time, int mesh);
void write_result(const PCGReturn &pcg_return, int mesh);

// PCG
PCGReturn pcg_solve(const LduMatrix &ldu_matrix, double * source, double *psi, int maxIter, double tolerance, double normfactor);
void ldu_to_csr(const LduMatrix &ldu_matrix, CsrMatrix &csr_matrix);
void csr_spmv(const CsrMatrix &csr_matrix, double *x, double *b);
void csr_precondition_spmv(const CsrMatrix &csr_matrix, double *vec, double *val, double *result);
void v_dot_product(const int nCells, const double *vec1, const double *vec2, double *result);
void v_sub_dot_product(const int nCells, const double *sub, const double *subed, const double *vec, double *result);
void pcg_init_precondition_csr (const CsrMatrix &csr_matrix, Precondition &pre);
void pcg_precondition_csr(const CsrMatrix &csr_matrix, const Precondition &pre, double *rAPtr, double *wAPtr);
double pcg_gsumMag(double *r, int size);
double pcg_gsumProd(double *z, double *r, int size);
void free_pcg(PCG &pcg);
void free_csr_matrix(CsrMatrix &csr_matrix);
void free_precondition(Precondition &pre);


#endif
