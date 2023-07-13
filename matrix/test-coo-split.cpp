#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>

#include "csr_matrix.cpp"
#include "matrix/coo_matrix.h"
#include "matrix/matrix_utils.h"
#include "pcg_def_opt.h"

int main(void) {
    srand(time(NULL));

    LduMatrix ldu_matrix;
    CsrMatrix csr_matrix;

    for (int _i = 0; _i < 25; ++_i) {
        int matrix_size = (_i + 1) * 64;
        fprintf(stderr, "matrix size: %dx%d\n", matrix_size, matrix_size);

        generate_ldu_matrix(ldu_matrix, matrix_size);

        ldu_to_csr(ldu_matrix, csr_matrix);
        auto splited = ldu_to_splited_coo(ldu_matrix, 64);

        free_ldu_matrix(ldu_matrix);

        auto plain_csr = plain_matrix_from_csr(csr_matrix);
        auto plain_splited = plain_matrix_from_splited_coo(splited);

        free_csr_matrix(csr_matrix);
        free_splited_coo(splited);

#ifdef SHOW_FIRST_100_ELEMENTS
        for (int i = 0; i < 100; ++i) {
            fprintf(stderr, "%f %f\n", plain_csr[i], plain_splited[i]);
        }
#endif

        assert(plain_csr == plain_splited);
    }

    return 0;
}
