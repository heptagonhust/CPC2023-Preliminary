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
    CscMatrix csc_matrix_from_ldu;
    CscMatrix csc_matrix_from_csr;

    for (int _i = 0; _i < 25; ++_i) {
        int matrix_size = (_i + 1) * 5;
        generate_ldu_matrix(ldu_matrix, matrix_size);
        ldu_to_csr(ldu_matrix, csr_matrix);
        ldu_to_csc(ldu_matrix, csc_matrix_from_ldu);
        csr_to_csc(csr_matrix, csc_matrix_from_csr);

        fprintf(
            stderr,
            "[INFO] matrix size = %d, csr size = %d, csc_from_ldu size = %d, compress_rate %.2f%%\n",
            matrix_size * matrix_size,
            csr_matrix.data_size,
            csc_matrix_from_ldu.data_size,
            csc_matrix_from_ldu.data_size * 100.0 / csr_matrix.data_size);

        auto plain_csc1 = plain_matrix_from_csc(csc_matrix_from_ldu);
        auto plain_csc2 = plain_matrix_from_csc(csc_matrix_from_csr);
        auto plain_csr = plain_matrix_from_csr(csr_matrix);

        free_csr_matrix(csr_matrix);
        free_csc_matrix(csc_matrix_from_csr);
        free_csc_matrix(csc_matrix_from_ldu);
        free_ldu_matrix(ldu_matrix);

        assert(
            plain_csc1 == plain_csc2 && plain_csc1 == plain_csr
            && plain_csr == plain_csc2);
    }

    return 0;
}
