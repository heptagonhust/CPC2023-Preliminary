
// cells : rows
// beta  : parameter
typedef struct {
	double *p_k;
	double *z_k1;
	double beta_k;
	int cells;
} MulAddPara;

typedef struct{
    double *m;
    double *r_k1;
		double *z_k1;
    int cells;
} MulPara;

typedef struct {
    double *r_k1;
    double *g;
    double *m;
    double *z_k1;
    int cells;
} SubMulPara;

typedef struct {
    double *r_k1;
    int cells;
    double normfactor;
    double tolerance;
		int *result;
		double *residual;
} ReducePara;

typedef struct {
    int cells;
    double *z_k1;
    double *r_k1;
    double *result;
} MulReduceZRPara;

typedef struct {
    int cells;
    double *p_k;
    double *Ax;
    double *result;
} MulReducepAxPara;

typedef struct {
    int cells;
    double *x;
    double *r;
    double *p;
    double *Ax;
    double alpha;
} UpdatexrPara;