#include "vector_slave.h"

#include <crts.h>
#include <slave.h>
#include <math.h>
#include <simd.h>
#include "pcg_def.h"
#include "vector_def.h"

// --------------------------------------------------------------------------------

inline void slave_MulAdd(double *p, double *z, double beta, int vec_num, int vec_start) {
    int front = 8 - vec_start % 8;
    int round = (vec_num - front) / 8;
    int end = (vec_num - front) % 8;
    int off = 0;
    doublev8 p8, z8, beta8;
    beta8 = simd_vcpyfd(beta);
    for (; off < front; ++off) {
        p[off] = z[off] + beta * p[off];
    }
    for (int i = 0; i < round; ++i, off += 8) {
        simd_load(p8, p + off);
        simd_load(z8, z + off);
        p8 = simd_vmad(beta8, p8, z8);
        simd_store(p8, p + off);
    }
    for (; off < vec_num; ++off) {
        p[off] = z[off] + beta * p[off];
    }
}

// --------------------------------------------------------------------------------

inline void slave_Mul(double *z, double *M_1, double *r, int vec_num) {
    double *aligned_z = (double *)(((uintptr_t)z + 63) & ~63);
    double *end_z = z + vec_num;
    double *round_z = (double *)((uintptr_t)end_z & ~63);
    for (; z < aligned_z; ++z, ++M_1, ++r) {
        *z = *M_1 * *r;
    }
    doublev8 z8, M_18, r8;
    for (; z < round_z; z += 8, M_1 += 8, r += 8) {
        simd_loadu(M_18, M_1);
        simd_loadu(r8, r);
        z8 = simd_vmuld(M_18, r8);
        simd_store(z8, z);
    }
    for (; z < end_z; ++z, ++M_1, ++r) {
        *z = *M_1 * *r;
    }
}

// --------------------------------------------------------------------------------

inline void slave_SubMul(double *z, double *r, double *g, double *M_1, int vec_num) {
    double *aligned_z = (double *)(((uintptr_t)z + 63) & ~63);
    double *end_z = z + vec_num;
    double *round_z = (double *)((uintptr_t)end_z & ~63);
    for (; z < aligned_z; ++z, ++r, ++g, ++M_1) {
        *z = (*r - *g) * *M_1;
    }
    doublev8 z8, r8, g8, M_18, tmp8;
    for (; z < round_z; z += 8, r += 8, g += 8, M_1 += 8) {
        simd_loadu(r8, r);
        simd_loadu(g8, g);
        simd_loadu(M_18, M_1);
        tmp8 = simd_vsubd(r8, g8);
        z8 = simd_vmuld(M_18, tmp8);
        simd_store(z8, z);
    }
    for (; z < end_z; ++z, ++r, ++g, ++M_1) {
        *z = (*r - *g) * *M_1;
    }
}

// --------------------------------------------------------------------------------

inline void slave_MulSub(double *g, double *z, double *M, int vec_num) {
    double *aligned_z = (double *)(((uintptr_t)z + 63) & ~63);
    double *end_z = z + vec_num;
    double *round_z = (double *)((uintptr_t)end_z & ~63);
    for (; z < aligned_z; ++g, ++z, ++M) {
        *g -= *z * *M;
    }
    doublev8 g8, z8, M8, tmp8;
    for (; z < round_z; g += 8, z += 8, M += 8) {
        simd_loadu(g8, g);
        simd_load(z8, z);
        simd_loadu(M8, M);
        tmp8 = simd_vnmad(z8, M8, g8);
        simd_storeu(tmp8, g);
    }
    for (; z < end_z; ++g, ++z, ++M) {
        *g -= *z * *M;
    }
}

// --------------------------------------------------------------------------------

inline double slave_Reduce(double *r, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    double *aligned_r = (double *)(((uintptr_t)r + 63) & ~63);
    double *end_r = r + vec_num;
    double *round_r = (double *)((uintptr_t)end_r & ~63);
    for (; r < aligned_r; ++r) {
        local_sum += fabs(*r);
    }

    union {
        double f64;
        int64_t i64;
    } mask;
    mask.i64 = 0x7fffffffffffffff;

    doublev8 local_sum8 = simd_vcpyfd(0.);
    union {
        int512 i64v8;
        doublev8 f64v8;
    } r8, mask8;
    mask8.f64v8 = simd_vcpyfd(mask.f64);
    for (; r < round_r; r += 8) {
        simd_load(r8.f64v8, r);
        r8.i64v8 &= mask8.i64v8;
        local_sum8 += r8.f64v8;
    }

    local_sum += simd_reduc_plusd(local_sum8);
    for (; r < end_r; ++r) {
        local_sum += fabs(*r);
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

inline double slave_MulReduceZR(double *r, double *z, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    doublev8 r8, z8, tmp, sum8 = simd_vcpyfd(0.);
    int round = vec_num / 8;
    for (int i = 0, off = 0; i < round; ++i, off += 8) {
        simd_loadu(r8, r + off);
        simd_loadu(z8, z + off);
        tmp = simd_vmad(r8, z8, sum8);
        sum8 = tmp;
    }
    local_sum = simd_reduc_plusd(sum8);
    for (int i = round * 8; i < vec_num; ++i) {
        local_sum += r[i] * z[i];
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

inline double slave_MulReducepAx(double *p, double *Ax, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    doublev8 p8, Ax8, tmp, local_sum8 = simd_vcpyfd(0.);
    int round = vec_num / 8;
    for (int i = 0, off = 0; i < round; ++i, off += 8) {
        simd_loadu(p8, p + off);
        simd_loadu(Ax8, Ax + off);
        tmp = simd_vmad(p8, Ax8, local_sum8);
        local_sum8 = tmp;
    }
    local_sum = simd_reduc_plusd(local_sum8);
    for (int i = round * 8; i < vec_num; ++i) {
        local_sum += p[i] * Ax[i];
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

inline void slave_Updatexr(double *x, double *r, double *p, double *Ax, double alpha, int vec_num) {
    double *aligned_p = (double *)(((uintptr_t)p + 63) & ~63);
    double *end_p = p + vec_num;
    double *round_p = (double *)((uintptr_t)end_p & ~63);
    for (; p < aligned_p; ++x, ++r, ++p, ++Ax) {
        *x += alpha * *p;
        *r -= alpha * *Ax;
    }
    doublev8 x8, r8, p8, Ax8, tmp8, alpha8;
    alpha8 = simd_vcpyfd(alpha);
    for (; p < round_p; x += 8, r += 8, p += 8, Ax += 8) {
        simd_loadu(x8, x);
        simd_loadu(r8, r);
        tmp8 = simd_vmad(alpha8, p8, x8);
        simd_storeu(tmp8, x);
        simd_load(p8, p);
        simd_loadu(Ax8, Ax);
        tmp8 = simd_vnmad(alpha8, Ax8, r8);
        simd_storeu(tmp8, r);
    }
    for (; p < end_p; ++x, ++r, ++p, ++Ax) {
        *x += alpha * *p;
        *r -= alpha * *Ax;
    }
}
