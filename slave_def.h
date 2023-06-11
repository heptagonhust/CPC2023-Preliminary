#ifndef _SLAVE_DEF_H_
#define _SLAVE_DEF_H_

#include <crts.h>

extern __thread_local crts_rply_t DMARply;
extern __thread_local unsigned int DMARplyCount;

/// 非阻塞 DMA 读
#define DMA_IGET(dst, src, len, rply) \
    CRTS_dma_iget(dst, src, len, rply); \
    DMARplyCount++

/// DMA 等待完成
#define DMA_WAIT CRTS_dma_wait_value(&DMARply, DMARplyCount)

/// 阻塞 DMA 读
#define DMA_GET(dst, src, len, rply) \
    CRTS_dma_iget(dst, src, len, rply); \
    DMARplyCount++; \
    DMA_WAIT

/// 阻塞 DMA 写
#define DMA_PUT(dst, src, len, rply) \
    CRTS_dma_iput(dst, src, len, rply); \
    DMARplyCount++; \
    DMA_WAIT

#endif