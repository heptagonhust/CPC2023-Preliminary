#pragma once

void v_dot_product(
    const int nCells,
    const double *vec1,
    const double *vec2,
    double *result) {
    for (int cell = 0; cell < nCells; cell++) {
        result[cell] = vec1[cell] * vec2[cell];
    }
}

void v_sub_dot_product(
    const int nCells,
    const double *sub,
    const double *subed,
    const double *vec,
    double *result) {
    for (int cell = 0; cell < nCells; cell++) {
        result[cell] = (sub[cell] - subed[cell]) * vec[cell];
    }
}
