#ifndef _SLAVE_DEF_H_
#define _SLAVE_DEF_H_

#include <crts.h>
#include "vector_def.h"

#define SLAVE_CORE_NUM 64

/// DMA 等待完成
#define DMA_WAIT(rply, cnt) CRTS_dma_wait_value(rply, cnt)

/// 非阻塞 DMA 读
#define DMA_IGET(dst, src, len, rply, cnt) \
    CRTS_dma_iget(dst, src, len, rply); \
    cnt++

/// 阻塞 DMA 读
#define DMA_GET(dst, src, len, rply, cnt) \
    CRTS_dma_iget(dst, src, len, rply); \
    cnt++; \
    DMA_WAIT(rply, cnt)

///非阻塞 DMA 写
#define DMA_IPUT(dst, src, len, rply, cnt) \
    CRTS_dma_iput(dst, src, len, rply); \
    cnt++;

#define DMA_IPUT_STRIDE(dst, src, len, bsize, stride, rply, cnt) \
    CRTS_dma_iput_stride(dst, src, len, bsize, stride, rply); \
    cnt++;

/// 阻塞 DMA 写
#define DMA_PUT(dst, src, len, rply, cnt) \
    CRTS_dma_iput(dst, src, len, rply); \
    cnt++; \
    DMA_WAIT(rply, cnt)


typedef struct {
    int maxIter;
    double normfactor;
    double tolerance;
    double *p;
    double *r;
    double *x;
    double *M;
    double *M_1;
    SpmvPara *spmv_para;
    Slave_task ntask[SLAVE_CORE_NUM];
    double *init_residual;
    double *final_residual;
    int *iter;
} MainLoopPara;


#endif