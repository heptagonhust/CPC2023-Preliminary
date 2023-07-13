#pragma once

#include <cstdlib>
#include <cstring>

#include "pcg_opt.h"
#include "pcg_def_opt.h"
#include "vector_utils.cpp"

void ldu_to_csr(const LduMatrix &ldu_matrix, CsrMatrix &csr_matrix) {
    csr_matrix.rows = ldu_matrix.cells;
    csr_matrix.data_size = 2 * ldu_matrix.faces + ldu_matrix.cells;
    csr_matrix.row_off = (int *)malloc((csr_matrix.rows + 1) * sizeof(int));
    csr_matrix.cols = (int *)malloc(csr_matrix.data_size * sizeof(int));
    csr_matrix.data = (double *)malloc(csr_matrix.data_size * sizeof(double));

    int row, col, offset;
    int *tmp = (int *)malloc((csr_matrix.rows + 1) * sizeof(int));

    csr_matrix.row_off[0] = 0;
    for (int i = 1; i < csr_matrix.rows + 1; i++)
        csr_matrix.row_off[i] = 1;

    for (int i = 0; i < ldu_matrix.faces; i++) {
        row = ldu_matrix.uPtr[i];
        col = ldu_matrix.lPtr[i];
        csr_matrix.row_off[row + 1]++;
        csr_matrix.row_off[col + 1]++;
    }

    for (int i = 0; i < ldu_matrix.cells; i++) {
        csr_matrix.row_off[i + 1] += csr_matrix.row_off[i];
    }

    memcpy(
        &tmp[0],
        &csr_matrix.row_off[0],
        (ldu_matrix.cells + 1) * sizeof(int));
    // lower
    for (int i = 0; i < ldu_matrix.faces; i++) {
        row = ldu_matrix.uPtr[i];
        col = ldu_matrix.lPtr[i];
        offset = tmp[row]++;
        csr_matrix.cols[offset] = col;
        csr_matrix.data[offset] = ldu_matrix.lower[i];
    }

    // diag
    for (int i = 0; i < ldu_matrix.cells; i++) {
        offset = tmp[i]++;
        csr_matrix.cols[offset] = i;
        csr_matrix.data[offset] = ldu_matrix.diag[i];
    }

    // upper
    for (int i = 0; i < ldu_matrix.faces; i++) {
        row = ldu_matrix.lPtr[i];
        col = ldu_matrix.uPtr[i];
        offset = tmp[row]++;
        csr_matrix.cols[offset] = col;
        csr_matrix.data[offset] = ldu_matrix.upper[i];
    }

    free(tmp);
}

// basic spmv, 需要负载均衡
void csr_spmv(const CsrMatrix &csr_matrix, double *vec, double *result) {
    for (int i = 0; i < csr_matrix.rows; i++) {
        int start = csr_matrix.row_off[i];
        int num = csr_matrix.row_off[i + 1] - csr_matrix.row_off[i];
        double temp = 0;
        for (int j = 0; j < num; j++) {
            temp +=
                vec[csr_matrix.cols[start + j]] * csr_matrix.data[start + j];
        }
        result[i] = temp;
    }
}

void csr_precondition_spmv(
    const CsrMatrix &csr_matrix,
    double *vec,
    double *val,
    double *result) {
    for (int i = 0; i < csr_matrix.rows; i++) {
        int start = csr_matrix.row_off[i];
        int num = csr_matrix.row_off[i + 1] - csr_matrix.row_off[i];
        double temp = 0;
        for (int j = 0; j < num; j++) {
            temp += vec[csr_matrix.cols[start + j]] * val[start + j];
        }
        result[i] = temp;
    }
}

// diagonal precondition, get matrix M^(-1) (diagonal matrix)
// pre_mat_val: 非对角元     : csr_matrix中元素
//              对角元素     : 0
// preD       : csr_matrix中对角元素的倒数
void pcg_init_precondition_csr(const CsrMatrix &csr_matrix, Precondition &pre) {
    for (int i = 0; i < csr_matrix.rows; i++) {
        for (int j = csr_matrix.row_off[i]; j < csr_matrix.row_off[i + 1];
             j++) {
            // get diagonal matrix
            if (csr_matrix.cols[j] == i) {
                pre.pre_mat_val[j] = 0.;
                pre.preD[i] = 1.0 / csr_matrix.data[j];
            } else {
                pre.pre_mat_val[j] = csr_matrix.data[j];
            }
        }
    }
}

// ? 存疑，循环用处?
void pcg_precondition_csr(
    const CsrMatrix &csr_matrix,
    const Precondition &pre,
    double *rAPtr,
    double *wAPtr) {
    double *gAPtr = (double *)malloc(csr_matrix.rows * sizeof(double));
    v_dot_product(csr_matrix.rows, pre.preD, rAPtr, wAPtr);
    memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
    for (int deg = 1; deg < 2; deg++) {
        // gAPtr = wAptr * pre.pre_mat_val; vec[rows] = matrix * vec[rows]
        csr_precondition_spmv(csr_matrix, wAPtr, pre.pre_mat_val, gAPtr);
        v_sub_dot_product(csr_matrix.rows, rAPtr, gAPtr, pre.preD, wAPtr);
        memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
    }
    free(gAPtr);
}

void free_csr_matrix(CsrMatrix &csr_matrix) {
    free(csr_matrix.cols);
    free(csr_matrix.data);
    free(csr_matrix.row_off);
}
