#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>

#include "csc_matrix.cpp"
#include "csr_matrix.cpp"
#include "matrix_utils.cpp"
#include "pcg_def.h"
#include "vector_utils.cpp"

int main(void) {
    srand(time(NULL));

    LduMatrix ldu_matrix;
    CsrMatrix csr_matrix;
    CscMatrix csc_matrix;

    for (int _i = 0; _i < 25; ++_i) {
        int matrix_size = (_i + 1) * 5;
        int vec_size = matrix_size;

        generate_ldu_matrix(ldu_matrix, matrix_size);
        ldu_to_csr(ldu_matrix, csr_matrix);
        ldu_to_csc(ldu_matrix, csc_matrix);
        free_ldu_matrix(ldu_matrix);

        double *vec = (double *)malloc(sizeof(double) * vec_size);
        double *csr_result = (double *)malloc(sizeof(double) * vec_size);
        double *csc_result = (double *)malloc(sizeof(double) * vec_size);
        for (int i = 0; i < vec_size; ++i) {
            vec[i] = (double)rand() / RAND_MAX * 20 - 10;
        }

        csc_spmv(csc_matrix, vec, csc_result);
        csr_spmv(csr_matrix, vec, csr_result);

        for (int i = 0; i < vec_size; ++i) {
            assert(csr_result[i] == csc_result[i]);
        }

        fprintf(stderr, "[INFO] matrix size %d\n", matrix_size);

        free(vec);
        free(csr_result);
        free(csc_result);
        free_csr_matrix(csr_matrix);
        free_csc_matrix(csc_matrix);
    }

    return 0;
}
