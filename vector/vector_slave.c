#include "vector_slave.h"

#include <crts.h>
#include <slave.h>
#include <math.h>
#include "pcg_def.h"
#include "vector_def.h"

// --------------------------------------------------------------------------------

void slave_MulAdd(double *p, double *z, double beta, int vec_num) {
    for (int i = 0; i < vec_num; i++) {
        p[i] = z[i] + beta * p[i];
    }
}

// --------------------------------------------------------------------------------

void slave_Mul(double *z, double *M_1, double *r, int vec_num) {
    for (int i = 0; i < vec_num; ++i) {
        z[i] = M_1[i] * r[i];
    }
}

// --------------------------------------------------------------------------------

void slave_SubMul(double *z, double *r, double *g, double *M_1, int vec_num) {
    for (int i = 0; i < vec_num; ++i) {
        z[i] = (r[i] - g[i]) * M_1[i];
    }
}

// --------------------------------------------------------------------------------

void slave_MulSub(double *g, double *z, double *M, int vec_num) {
    for (int i = 0; i < vec_num; ++i) {
        g[i] = g[i] - z[i] * M[i];
    }
}

// --------------------------------------------------------------------------------

double slave_Reduce(double *r, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    for (int i = 0; i < vec_num; ++i) {
        local_sum += fabs(r[i]);
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

double slave_MulReduceZR(double *r, double *z, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    for (int i = 0; i < vec_num; ++i) {
        local_sum += r[i] * z[i];
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

double slave_MulReducepAx(double *p, double *Ax, int vec_num, double *reducebuf) {
    double local_sum = 0, sum = 0;
    for (int i = 0; i < vec_num; ++i) {
        local_sum += p[i] * Ax[i];
    }
    CRTS_scoll_redurt(&local_sum, &sum, 1, CRTS_double, OP_add, reducebuf, 64);
    return sum;
}

// --------------------------------------------------------------------------------

void slave_Updatexr(double *x, double *r, double *p, double *Ax, double alpha, int vec_num) {
    for (int i = 0; i < vec_num; ++i) {
        x[i] = x[i] + alpha * p[i];
        r[i] = r[i] - alpha * Ax[i];
    }
}
