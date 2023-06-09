#pragma once

#include <cstdlib>
#include <cstring>

#include "matrix_utils.cpp"
#include "pcg.h"
#include "pcg_def.h"

void ldu_to_csc(const LduMatrix &ldu_matrix, CscMatrix &csc_matrix) {
    csc_matrix.cols = ldu_matrix.cells;
    csc_matrix.data_size = count_ldu_matrix_nonzero_elements(ldu_matrix);
    csc_matrix.col_off = (int *)malloc((csc_matrix.cols + 1) * sizeof(int));
    csc_matrix.rows = (int *)malloc(csc_matrix.data_size * sizeof(int));
    csc_matrix.data = (double *)malloc(csc_matrix.data_size * sizeof(double));

    int row, col, offset;
    int *tmp = (int *)malloc((csc_matrix.cols + 1) * sizeof(int));

    csc_matrix.col_off[0] = 0;
    for (int i = 1; i < ldu_matrix.cells + 1; i++) {
        if (ldu_matrix.diag[i - 1] != 0.0) {
            csc_matrix.col_off[i] = 1;
        } else {
            csc_matrix.col_off[i] = 0;
        }
    }

    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.lower[i] != 0.0) {
            col = ldu_matrix.lPtr[i];
            csc_matrix.col_off[col + 1]++;
        }
    }

    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.upper[i] != 0.0) {
            row = ldu_matrix.uPtr[i];
            csc_matrix.col_off[row + 1]++;
        }
    }

    for (int i = 0; i < ldu_matrix.cells; i++) {
        csc_matrix.col_off[i + 1] += csc_matrix.col_off[i];
    }

    memcpy(
        &tmp[0],
        &csc_matrix.col_off[0],
        (ldu_matrix.cells + 1) * sizeof(int));

    // lower
    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.lower[i] != 0.0) {
            row = ldu_matrix.uPtr[i];
            col = ldu_matrix.lPtr[i];
            offset = tmp[col]++;
            csc_matrix.rows[offset] = row;
            csc_matrix.data[offset] = ldu_matrix.lower[i];
        }
    }

    // diag
    for (int i = 0; i < ldu_matrix.cells; i++) {
        if (ldu_matrix.diag[i] != 0.0) {
            offset = tmp[i]++;
            csc_matrix.rows[offset] = i;
            csc_matrix.data[offset] = ldu_matrix.diag[i];
        }
    }

    // upper
    for (int i = 0; i < ldu_matrix.faces; i++) {
        if (ldu_matrix.upper[i] != 0.0) {
            row = ldu_matrix.lPtr[i];
            col = ldu_matrix.uPtr[i];
            offset = tmp[col]++;
            csc_matrix.rows[offset] = row;
            csc_matrix.data[offset] = ldu_matrix.upper[i];
        }
    }

    free(tmp);
}

// // basic spmv, 需要负载均衡
// void csc_spmv(const CscMatrix &csc_matrix, double *vec, double *result) {
//     for (int i = 0; i < csr_matrix.rows; i++) {
//         int start = csr_matrix.row_off[i];
//         int num = csr_matrix.row_off[i + 1] - csr_matrix.row_off[i];
//         double temp = 0;
//         for (int j = 0; j < num; j++) {
//             temp +=
//                 vec[csr_matrix.cols[start + j]] * csr_matrix.data[start + j];
//         }
//         result[i] = temp;
//     }
// }
//
// void csc_precondition_spmv(
//     const CscMatrix &csc_matrix,
//     double *vec,
//     double *val,
//     double *result) {
//     for (int i = 0; i < csr_matrix.rows; i++) {
//         int start = csr_matrix.row_off[i];
//         int num = csr_matrix.row_off[i + 1] - csr_matrix.row_off[i];
//         double temp = 0;
//         for (int j = 0; j < num; j++) {
//             temp += vec[csr_matrix.cols[start + j]] * val[start + j];
//         }
//         result[i] = temp;
//     }
// }
//
// // diagonal precondition, get matrix M^(-1) (diagonal matrix)
// // pre_mat_val: 非对角元     : csr_matrix中元素
// //              对角元素     : 0
// // preD       : csr_matrix中对角元素的倒数
// void pcg_init_precondition_csc(const CscMatrix &csc_matrix, Precondition &pre) {
//     for (int i = 0; i < csr_matrix.rows; i++) {
//         for (int j = csr_matrix.row_off[i]; j < csr_matrix.row_off[i + 1];
//              j++) {
//             // get diagonal matrix
//             if (csr_matrix.cols[j] == i) {
//                 pre.pre_mat_val[j] = 0.;
//                 pre.preD[i] = 1.0 / csr_matrix.data[j];
//             } else {
//                 pre.pre_mat_val[j] = csr_matrix.data[j];
//             }
//         }
//     }
// }
//
// // ? 存疑，循环用处?
// void pcg_precondition_csc(
//     const CscMatrix &csc_matrix,
//     const Precondition &pre,
//     double *rAPtr,
//     double *wAPtr) {
//     double *gAPtr = (double *)malloc(csr_matrix.rows * sizeof(double));
//     v_dot_product(csr_matrix.rows, pre.preD, rAPtr, wAPtr);
//     memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
//     for (int deg = 1; deg < 2; deg++) {
//         // gAPtr = wAptr * pre.pre_mat_val; vec[rows] = matrix * vec[rows]
//         csr_precondition_spmv(csr_matrix, wAPtr, pre.pre_mat_val, gAPtr);
//         v_sub_dot_product(csr_matrix.rows, rAPtr, gAPtr, pre.preD, wAPtr);
//         memset(gAPtr, 0, csr_matrix.rows * sizeof(double));
//     }
//     free(gAPtr);
// }

void free_csc_matrix(CscMatrix &mtx) {
    free(mtx.rows);
    free(mtx.data);
    free(mtx.col_off);
}
