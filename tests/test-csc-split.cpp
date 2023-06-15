#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>

#include "csc_matrix.cpp"
#include "matrix_utils.cpp"
#include "pcg_def.h"

int main(void) {
    srand(time(NULL));

    LduMatrix ldu_matrix;
    CscMatrix csc_matrix;

    for (int _i = 0; _i < 25; ++_i) {
        int matrix_size = (_i + 1) * 64;
        fprintf(stderr, "matrix size: %dx%d\n", matrix_size, matrix_size);
        generate_ldu_matrix(ldu_matrix, matrix_size);
        ldu_to_csc(ldu_matrix, csc_matrix);
        free_ldu_matrix(ldu_matrix);

        auto splited = split_csc_matrix(csc_matrix, 64);
        for (int i = 0; i < splited.chunk_num; ++i) {
            csc_chunk_pack(*splited.chunks[i]);
            csc_chunk_unpack(*splited.chunks[i]);
        }

        auto plain_csc = plain_matrix_from_csc(csc_matrix);
        auto plain_splited = plain_matrix_from_splited_csc_matrix(splited);

        free_csc_matrix(csc_matrix);
        free_packed_splited_csc_matrix(&splited);

        assert(plain_csc == plain_splited);
    }

    return 0;
}
