#ifndef _PCG_DEF_H_
#define _PCG_DEF_H_


struct LduMatrix {
	double *upper;
	double *lower;
	double *diag;
	int *uPtr;
	int *lPtr;
	int faces;
	int cells;
};

struct CsrMatrix {
    int rows;
    int *row_off;
    int *cols;
    double *data;
    int data_size;
};

struct PCG {
    double *r;
    double *z;
    double *p;
    double *Ax;
    double sumprod;
    double sumprod_old;
    double residual;
    double alpha;
    double beta;

    double *x;
    double *source;
};

struct Precondition{
    double *pre_mat_val;
    double *preD;
};

struct PCGReturn{
    double residual;
    int iter;
};

#include <time.h>
#include <stdio.h>
#define INFO(M, ...) {  time_t t; \
                        struct tm *tmif; \
                        t = time(NULL); \
                        tmif = localtime(&t); \
                        printf("[%d-%2d-%2d] [%2d:%2d:%2d] [INFO] " M "",  \
                                tmif->tm_year + 1900, \
                                tmif->tm_mon + 1, \
                                tmif->tm_mday, \
                                tmif->tm_hour, \
                                tmif->tm_min, \
                                tmif->tm_sec, ##__VA_ARGS__); \
                    }

#endif
