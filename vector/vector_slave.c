#include "vector_slave.h"

#include <crts.h>
#include <slave.h>
#include <math.h>
#include <simd.h>
#include "pcg_def.h"
#include "vector_def.h"

#define SIMD_WIDTH 8
__thread_local double simd_arr[SIMD_WIDTH] __attribute__((aligned(64)));

// --------------------------------------------------------------------------------

inline void slave_MulAdd(double *p, double *z, double beta, int vec_num, int vec_start) {
    int front = 8 - vec_start % 8;
    int round = (vec_num - front) / 8;
    int end = (vec_num - front) % 8;
    int off = 0;
    doublev8 p8, z8, beta8;
    for (int i = 0; i < 8; ++i) simd_arr[i] = beta;
    simd_load(beta8, simd_arr);
    for (; off < front; ++off) {
        p[off] = z[off] + beta * p[off];
    }
    for (int i = 0; i < round; ++i, off += 8) {
        simd_load(p8, p + off);
        simd_load(z8, z + off);
        p8 = z8 + beta8 * p8;
        simd_store(p8, p + off);
    }
    for (; off < vec_num; ++off) {
        p[off] = z[off] + beta * p[off];
    }
}

// --------------------------------------------------------------------------------

inline void slave_Mul(double *z, double *M_1, double *r, int vec_num) {
    for (int i = 0; i < vec_num; ++i) {
        z[i] = M_1[i] * r[i];
    }
}

// --------------------------------------------------------------------------------

inline void slave_SubMul(double *z, double *r, double *g, double *M_1, int vec_num) {
    for (int i = 0; i < vec_num; ++i) {
        z[i] = (r[i] - g[i]) * M_1[i];
    }
}

// --------------------------------------------------------------------------------

inline void slave_MulSub(double *g, double *z, double *M, int vec_num) {
    for (int i = 0; i < vec_num; ++i) {
        g[i] = g[i] - z[i] * M[i];
    }
}

// --------------------------------------------------------------------------------

inline double slave_Reduce(double *r, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    for (int i = 0; i < vec_num; ++i) {
        local_sum += fabs(r[i]);
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

inline double slave_MulReduceZR(double *r, double *z, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    for (int i = 0; i < vec_num; ++i) {
        local_sum += r[i] * z[i];
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

inline double slave_MulReducepAx(double *p, double *Ax, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    for (int i = 0; i < vec_num; ++i) {
        local_sum += p[i] * Ax[i];
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

inline void slave_Updatexr(double *x, double *r, double *p, double *Ax, double alpha, int vec_num) {
    for (int i = 0; i < vec_num; ++i) {
        x[i] = x[i] + alpha * p[i];
        r[i] = r[i] - alpha * Ax[i];
    }
}
