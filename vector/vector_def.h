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
    void slave_MulAdd(MulAddPara *para);

    //! result[i] = vec0[i] * vec1[i]
    void slave_Mul(MulPara *para);

    //! result[i] = (vec0[i] - vec1[i]) * vec2[i]
    void slave_SubMul(SubMulPara *para);

    //! vec0[i] = vec0[i] - vec1[i] * vec2[i]
    void slave_MulSub(MulSubPara *para);

    //! result += abs(vec[i])
    void slave_Reduce(ReducePara *para);

    //! result += vec0[i] * vec1[i]
    void slave_MulReduceZR(MulReduceZRPara *para);

    //! result += vec2[i] * vec3[i]
    void slave_MulReducepAx(MulReducepAxPara *para);

    void slave_Updatexr(UpdatexrPara *para);

#ifdef __cplusplus
}
#endif
