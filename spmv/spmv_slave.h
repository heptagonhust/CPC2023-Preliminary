#ifndef _SPMV_H_
#define _SPMV_H_

#include "spmv_def.h"
#include <crts.h>
#include <slave.h>
#include "slave_def.h"

typedef struct {
    void *buff[2];
    int size[2];
    int in_use;
} DoubleBuffering;

extern void *slave_double_buffering_get(DoubleBuffering *buff, int size);
extern void slave_double_buffering_free(DoubleBuffering *buff);

extern void slave_csc_spmv(void *para);

#endif