#ifndef _SPMV_H_
#define _SPMV_H_

#include "spmv_def.h"
#include "slave_def.h"
#include <crts.h>
#include <slave.h>
#include <string.h>

typedef struct {
    double *buff[2];
    int size;
    int in_use;
} DoubleBuffering;


#ifdef __cplusplus
extern "C" {
#endif

  void slave_coo_spmv(CooChunk *chunk, double *vec, double **vec_list, double *result, DoubleBuffering *buff, uint16_t non_0_block_num, uint16_t *non_0_block_idx);

#ifdef __cplusplus
}
#endif


void slave_double_buffering_new(DoubleBuffering *buff, int size);
void* slave_double_buffering_get(DoubleBuffering *buff);
void slave_double_buffering_free(DoubleBuffering *buff);

#endif