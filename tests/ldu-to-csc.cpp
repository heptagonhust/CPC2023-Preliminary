#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <ctime>
#include <vector>

#include "../pcg_def.h"
#include "../pcg.h"

void csr_to_csc(const CsrMatrix &csr_matrix, CscMatrix &csc_matrix) {
    csc_matrix.cols = csr_matrix.rows;
    csc_matrix.data_size = csr_matrix.data_size;
    csc_matrix.col_off = (int *)malloc((csc_matrix.cols + 1)*sizeof(int));
    csc_matrix.rows = (int *)malloc(csc_matrix.data_size*sizeof(int));
    csc_matrix.data = (double *)malloc(csc_matrix.data_size*sizeof(double));

    for (int i = 0; i < csc_matrix.cols + 1; i++) {
        csc_matrix.col_off[i] = 0;
    }

    for (int i = 0; i < csr_matrix.data_size; i++) {
        csc_matrix.col_off[csr_matrix.cols[i] + 1]++;
    }

    int last_column_size = csc_matrix.col_off[csc_matrix.cols -1];

    for (int i = 0; i < csc_matrix.cols; i++) {
        csc_matrix.col_off[i + 1] += csc_matrix.col_off[i];
    }

    for (int i = 0; i < csr_matrix.rows; i++) {
        for (int j = csr_matrix.row_off[i]; j < csr_matrix.row_off[i + 1]; j++) {
            int col = csr_matrix.cols[j];
            int dest = csc_matrix.col_off[col];

            csc_matrix.rows[dest] = i;
            csc_matrix.data[dest] = csr_matrix.data[j];

            csc_matrix.col_off[col]++;
        }
    }
    for (int i = 0; i < csc_matrix.cols; i++) {
        csc_matrix.col_off[i] -= csc_matrix.col_off[i + 1] - csc_matrix.col_off[i];
    }
    csc_matrix.col_off[csc_matrix.cols - 1] -= last_column_size;
}

void ldu_to_csc(const LduMatrix &ldu_matrix, CscMatrix &csc_matrix) {
    csc_matrix.cols = ldu_matrix.cells;
    csc_matrix.data_size = 2*ldu_matrix.faces + ldu_matrix.cells;
    csc_matrix.col_off = (int *)malloc((csc_matrix.cols + 1)*sizeof(int));
    csc_matrix.rows = (int *)malloc(csc_matrix.data_size*sizeof(int));
    csc_matrix.data = (double *)malloc(csc_matrix.data_size*sizeof(double));

    int row, col, offset;
    int *tmp = (int *)malloc((csc_matrix.cols + 1)*sizeof(int));

    csc_matrix.col_off[0] = 0;
    for(int i = 1; i < csc_matrix.cols + 1; i++)
        csc_matrix.col_off[i] = 1;

    for(int i = 0; i < ldu_matrix.faces; i++){
        row     = ldu_matrix.uPtr[i] ;
        col = ldu_matrix.lPtr[i] ;
        csc_matrix.col_off[row+1]++;
        csc_matrix.col_off[col+1]++;
    }

    for(int i = 0;i< ldu_matrix.cells; i++){
        csc_matrix.col_off[i+1] += csc_matrix.col_off[i];
    }

    memcpy(&tmp[0], &csc_matrix.col_off[0], (ldu_matrix.cells + 1)*sizeof(int));
    // lower
    for(int i = 0; i < ldu_matrix.faces; i++ ){
        row = ldu_matrix.uPtr[i];
        col = ldu_matrix.lPtr[i];
        offset = tmp[col]++;
        csc_matrix.rows[offset] = row;
        csc_matrix.data[offset] = ldu_matrix.lower[i];
    }

    // diag
    for(int i = 0; i < ldu_matrix.cells; i++){
        offset = tmp[i]++;
        csc_matrix.rows[offset] = i;
        csc_matrix.data[offset] = ldu_matrix.diag[i];
    }

    // upper
    for(int i = 0; i < ldu_matrix.faces; i++){
        row = ldu_matrix.lPtr[i];
        col = ldu_matrix.uPtr[i];
        offset = tmp[col]++;
        csc_matrix.rows[offset] = row;
        csc_matrix.data[offset] = ldu_matrix.upper[i];
    }

    free(tmp);
}

void ldu_to_csr(const LduMatrix &ldu_matrix, CsrMatrix &csr_matrix) {
    csr_matrix.rows = ldu_matrix.cells;
    csr_matrix.data_size = 2*ldu_matrix.faces + ldu_matrix.cells;
    csr_matrix.row_off = (int *)malloc((csr_matrix.rows + 1)*sizeof(int));
    csr_matrix.cols = (int *)malloc(csr_matrix.data_size*sizeof(int));
    csr_matrix.data = (double *)malloc(csr_matrix.data_size*sizeof(double));

    int row, col, offset;
    int *tmp = (int *)malloc((csr_matrix.rows + 1)*sizeof(int));

    csr_matrix.row_off[0] = 0;
    for(int i = 1; i < csr_matrix.rows + 1; i++)
        csr_matrix.row_off[i] = 1;

    for(int i = 0; i < ldu_matrix.faces; i++){
        row	= ldu_matrix.uPtr[i] ;
        col = ldu_matrix.lPtr[i] ;
        csr_matrix.row_off[row+1]++;
        csr_matrix.row_off[col+1]++;
    }

    for(int i = 0;i< ldu_matrix.cells; i++){
        csr_matrix.row_off[i+1] += csr_matrix.row_off[i];
    }

    memcpy(&tmp[0], &csr_matrix.row_off[0], (ldu_matrix.cells + 1)*sizeof(int));
    // lower
    for(int i = 0; i < ldu_matrix.faces; i++ ){
        row = ldu_matrix.uPtr[i];
        col = ldu_matrix.lPtr[i];
        offset = tmp[row]++;
        csr_matrix.cols[offset] = col;
        csr_matrix.data[offset] = ldu_matrix.lower[i];
    }

    // diag
    for(int i = 0; i < ldu_matrix.cells; i++){
        offset = tmp[i]++;
        csr_matrix.cols[offset] = i;
        csr_matrix.data[offset] = ldu_matrix.diag[i];
    }

    // upper
    for(int i = 0; i < ldu_matrix.faces; i++){
        row = ldu_matrix.lPtr[i];
        col = ldu_matrix.uPtr[i];
        offset = tmp[row]++;
        csr_matrix.cols[offset] = col;
        csr_matrix.data[offset] = ldu_matrix.upper[i];
    }

    free(tmp);
}

double rand_kernel(void) {
    if (rand() & 1) {
        return 0.;
    }

    return (double)rand() / RAND_MAX * 10;
}

void generate_ldu_matrix(LduMatrix &ldu_matrix, int size) {
    ldu_matrix.cells = size;
    ldu_matrix.faces = size * (size - 1) / 2;
    ldu_matrix.lPtr = (int *)malloc(ldu_matrix.faces * sizeof(int));
    ldu_matrix.uPtr = (int *)malloc(ldu_matrix.faces * sizeof(int));
    ldu_matrix.lower = (double *)malloc(ldu_matrix.faces * sizeof(double));
    ldu_matrix.upper = (double *)malloc(ldu_matrix.faces * sizeof(double));
    ldu_matrix.diag = (double *)malloc(ldu_matrix.cells * sizeof(double));

    int k = 0;
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            ldu_matrix.lPtr[k] = i;
            ldu_matrix.uPtr[k] = j;
            ldu_matrix.lower[k] = rand_kernel();
            ldu_matrix.upper[k] = rand_kernel();
            k++;
        }
        ldu_matrix.diag[i] = rand_kernel();
    }
}

std::vector<double> plain_matrix_from_csr(const CsrMatrix &csr_matrix) {
    std::vector<double> result;

    for (int i = 0; i < csr_matrix.rows; i++) {
        int k = csr_matrix.row_off[i];
        for (int j = 0; j < csr_matrix.rows; j++) {
            if (k < csr_matrix.row_off[i + 1] && csr_matrix.cols[k] == j) {
                result.push_back(csr_matrix.data[k]);
                k++;
            } else {
                result.push_back(0.);
            }
        }
    }

    return result;
}

std::vector<double> plain_matrix_from_csc(const CscMatrix &csc_matrix) {
    std::vector<double> result;

    for (int i = 0; i < csc_matrix.cols; i++) {
        for (int j = 0; j < csc_matrix.cols; j++) {
            int k = csc_matrix.col_off[j];
            bool found = false;
            while (k < csc_matrix.col_off[j + 1]) {
                if (csc_matrix.rows[k] == i) {
                    result.push_back(csc_matrix.data[k]);
                    found = true;
                    break;
                }
                k++;
            }
            if (!found) {
                result.push_back(0.);
            }
        }
    }

    return result;
}

bool csc_check(const CscMatrix &lhs, const CscMatrix &rhs) {
#define DO_CHECK(a, b) do { \
    if ((a) != (b)) { \
        fprintf(stderr, "[ERR] %s != %s\n", #a, #b); \
        return false; \
    } \
} while (0)

    DO_CHECK(lhs.cols, rhs.cols);
    DO_CHECK(lhs.data_size, rhs.data_size);

    for (int i = 0; i < lhs.cols + 1; ++i) {
        DO_CHECK(lhs.col_off[i], rhs.col_off[i]);
    }

    for (int i = 0; i < lhs.data_size; ++i) {
        DO_CHECK(lhs.rows[i], rhs.rows[i]);
        DO_CHECK(lhs.data[i], rhs.data[i]);
    }

    return true;
}

int main(void) {
    srand(time(NULL));

    LduMatrix ldu_matrix;
    CsrMatrix csr_matrix;
    CscMatrix csc_matrix_from_ldu;
    CscMatrix csc_matrix_from_csr;
    generate_ldu_matrix(ldu_matrix, 5);
    ldu_to_csr(ldu_matrix, csr_matrix);
    ldu_to_csc(ldu_matrix, csc_matrix_from_ldu);
    csr_to_csc(csr_matrix, csc_matrix_from_csr);

    auto plain_csc1 = plain_matrix_from_csc(csc_matrix_from_ldu);
    auto plain_csc2 = plain_matrix_from_csc(csc_matrix_from_csr);
    auto plain_csr = plain_matrix_from_csr(csr_matrix);

    assert(plain_csc1 == plain_csc2 && plain_csc1 == plain_csr && plain_csr == plain_csc2);

    return 0;
}
