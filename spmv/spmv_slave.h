#ifndef _SPMV_H_
#define _SPMV_H_

#include "spmv_def.h"
#include "slave_def.h"
#include <crts.h>
#include <slave.h>
#include <string.h>

typedef struct {
    void *buff[2];
    int size;
    int in_use;
} DoubleBuffering;

void slave_csc_spmv(SpmvPara *para_mp);

#endif