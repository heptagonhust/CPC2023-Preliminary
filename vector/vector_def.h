#pragma once

typedef struct {
    int col_start;
    int col_num;
} Slave_task;

typedef struct {
    double *p_k;
    double *z_k1;
    double beta_k;
    int cells;
    Slave_task task[64];
} MulAddPara;

typedef struct {
    double *m;
    double *r_k1;
    double *z_k1;
    int cells;
    Slave_task task[64];
} MulPara;

typedef struct {
    double *r_k1;
    double *g;
    double *m;
    double *z_k1;
    int cells;
    Slave_task task[64];
} SubMulPara;

typedef struct {
    double *g;
    double *z_k1;
    double *m;
    int cells;
    Slave_task task[64];
} MulSubPara;

typedef struct {
    double *r_k1;
    int cells;
    double normfactor;
    double tolerance;
    int *result;
    double *residual;
    Slave_task task[64];
} ReducePara;

typedef struct {
    int cells;
    double *z_k1;
    double *r_k1;
    double *result;
    Slave_task task[64];
} MulReduceZRPara;

typedef struct {
    int cells;
    double *p_k;
    double *Ax;
    double *result;
    Slave_task task[64];
} MulReducepAxPara;

typedef struct {
    int cells;
    double *x;
    double *r;
    double *p;
    double *Ax;
    double alpha;
    Slave_task task[64];
} UpdatexrPara;

#ifdef __cplusplus
extern "C" {
#endif
    //! vec1 = vec0 + scalar * vec1
    void slave_MulAdd(double *p, double *z, double beta, int vec_num, int vec_start);

    //! result[i] = vec0[i] * vec1[i]
    void slave_Mul(double *z, double *M_1, double *r, int vec_num);

    //! result[i] = (vec0[i] - vec1[i]) * vec2[i]
    void slave_SubMul(double *z, double *r, double *g, double *M_1, int vec_num);

    //! vec0[i] = vec0[i] - vec1[i] * vec2[i]
    void slave_MulSub(double *g, double *z, double *M, int vec_num);

    //! result += abs(vec[i])
    double slave_Reduce(double *r, int vec_num, double *reducebuf);

    //! result += vec0[i] * vec1[i]
    double slave_MulReduceZR(double *r, double *z, int vec_num, double *reducebuf);

    //! result += vec2[i] * vec3[i]
    double slave_MulReducepAx(double *p, double *Ax, int vec_num, double *reducebuf);

    void slave_Updatexr(double *x, double *r, double *p, double *Ax, double alpha, int vec_num);

#ifdef __cplusplus
}
#endif
