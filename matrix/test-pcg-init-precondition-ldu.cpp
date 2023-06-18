#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>

#include "csr_matrix.cpp"
#include "matrix/coo_matrix.h"
#include "matrix/matrix_utils.h"
#include "pcg_def.h"

void pcg_init_precondition_csr(
    const CsrMatrix &csr_matrix,
    double *pre,
    double *M) {
    for (int i = 0; i < csr_matrix.rows; i++) {
        for (int j = csr_matrix.row_off[i]; j < csr_matrix.row_off[i + 1];
             j++) {
            if (csr_matrix.cols[j] == i) {
                pre[i] = 1.0 / csr_matrix.data[j];
                M[i] = csr_matrix.data[j];
            }
        }
    }
}

void pcg_init_precondition_ldu(
    const LduMatrix &ldu_matrix,
    double *pre,
    double *M) {
    for (int i = 0; i < ldu_matrix.cells; i++) {
        pre[i] = 1.0 / ldu_matrix.diag[i];
        M[i] = ldu_matrix.diag[i];
    }
}

int main(void) {
    srand(time(NULL));

    LduMatrix ldu_matrix;
    CsrMatrix csr_matrix;

    for (int _i = 0; _i < 25; ++_i) {
        int matrix_size = (_i + 1) * 64;
        fprintf(stderr, "matrix size: %dx%d\n", matrix_size, matrix_size);

        generate_ldu_matrix(ldu_matrix, matrix_size);

        ldu_to_csr(ldu_matrix, csr_matrix);
        double *result_pre_csr =
            (double *)malloc(sizeof(double) * ldu_matrix.cells);
        double *result_M_csr =
            (double *)malloc(sizeof(double) * ldu_matrix.cells);
        double *result_pre_ldu =
            (double *)malloc(sizeof(double) * ldu_matrix.cells);
        double *result_M_ldu =
            (double *)malloc(sizeof(double) * ldu_matrix.cells);

        pcg_init_precondition_csr(csr_matrix, result_pre_csr, result_M_csr);
        pcg_init_precondition_ldu(ldu_matrix, result_pre_ldu, result_M_ldu);

        for (int i = 0; i < csr_matrix.rows; ++i) {
            assert(result_pre_csr[i] == result_pre_ldu[i]);
            assert(result_M_csr[i] == result_M_ldu[i]);
        }

        free_ldu_matrix(ldu_matrix);
        free_csr_matrix(csr_matrix);
        free(result_pre_csr);
        free(result_M_csr);
        free(result_pre_ldu);
        free(result_M_ldu);
    }

    return 0;
}
