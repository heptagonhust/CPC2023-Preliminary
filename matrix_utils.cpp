#pragma once

#include <vector>
#include "pcg_def.h"

int count_ldu_matrix_nonzero_elements(const LduMatrix &ldu_matrix) {
    int result = 0;
    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.upper[i] != 0.0) {
            result += 1;
        }
    }

    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.lower[i] != 0.0) {
            result += 1;
        }
    }

    for (int i = 0; i < ldu_matrix.cells; i++) {
        if (ldu_matrix.diag[i] != 0.0) {
            result += 1;
        }
    }

    return result;
}

void csr_to_csc(const CsrMatrix &csr_matrix, CscMatrix &csc_matrix) {
    csc_matrix.cols = csr_matrix.rows;
    csc_matrix.data_size = csr_matrix.data_size;
    csc_matrix.col_off = (int *)malloc((csc_matrix.cols + 1) * sizeof(int));
    csc_matrix.rows = (int *)malloc(csc_matrix.data_size * sizeof(int));
    csc_matrix.data = (double *)malloc(csc_matrix.data_size * sizeof(double));

    for (int i = 0; i < csc_matrix.cols + 1; i++) {
        csc_matrix.col_off[i] = 0;
    }

    for (int i = 0; i < csr_matrix.data_size; i++) {
        csc_matrix.col_off[csr_matrix.cols[i] + 1]++;
    }

    int last_column_size = csc_matrix.col_off[csc_matrix.cols - 1];

    for (int i = 0; i < csc_matrix.cols; i++) {
        csc_matrix.col_off[i + 1] += csc_matrix.col_off[i];
    }

    for (int i = 0; i < csr_matrix.rows; i++) {
        for (int j = csr_matrix.row_off[i]; j < csr_matrix.row_off[i + 1];
             j++) {
            int col = csr_matrix.cols[j];
            int dest = csc_matrix.col_off[col];

            csc_matrix.rows[dest] = i;
            csc_matrix.data[dest] = csr_matrix.data[j];

            csc_matrix.col_off[col]++;
        }
    }
    for (int i = 0; i < csc_matrix.cols; i++) {
        csc_matrix.col_off[i] -=
            csc_matrix.col_off[i + 1] - csc_matrix.col_off[i];
    }
    csc_matrix.col_off[csc_matrix.cols - 1] -= last_column_size;
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

void free_ldu_matrix(LduMatrix &mtx) {
    free(mtx.upper);
    free(mtx.lower);
    free(mtx.diag);
    free(mtx.lPtr);
    free(mtx.uPtr);
}
