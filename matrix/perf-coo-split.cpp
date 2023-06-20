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

int main(void) {
    srand(time(NULL));
    LduMatrix ldu_matrix;
    int matrix_size = 25 * 64;
    generate_ldu_matrix(ldu_matrix, matrix_size);
    auto splited = ldu_to_splited_coo(ldu_matrix, 64);
    free_ldu_matrix(ldu_matrix);
    return 0;
}
