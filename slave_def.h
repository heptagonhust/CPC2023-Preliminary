
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