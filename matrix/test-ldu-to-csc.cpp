#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>

#include "matrix/csc_matrix.h"
#include "csr_matrix.cpp"
#include "matrix/matrix_utils.h"
#include "pcg_def_opt.h"
#include "vector_utils.cpp"

bool check_csc_matrix_row_number_order(CscMatrix &mtx) {
    for (int i = 0; i < mtx.cols; ++i) {
        int size = mtx.col_off[i + 1] - mtx.col_off[i];
        int *rows = mtx.rows + mtx.col_off[i];
        for (int j = 1; j < size; ++j) {
            if (rows[j - 1] >= rows[j]) {
                return false;
            }
        }
    }
    return true;
}
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

        assert(check_csc_matrix_row_number_order(csc_matrix_from_csr));

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
