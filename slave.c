#include <slave.h>
#include "pcg_def.h"
#include "spmv.h"

#include <crts.h>

typedef struct{
	double *p;
	double *z;
	double beta;
	int cells;
} Para;

typedef struct {
    // TODO
} CscChunkPara;

#define dataBufferSize 2000
__thread_local crts_rply_t DMARply = 0;
__thread_local unsigned int DMARplyCount = 0;
__thread_local double p[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double z[dataBufferSize] __attribute__ ((aligned(64)));

void slave_example(Para* para){
	Para slavePara;
	//接收结构体数据
	CRTS_dma_iget(&slavePara, para, sizeof(Para), &DMARply);
	DMARplyCount++;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
	double beta = slavePara.beta;
	int cells = slavePara.cells;
	
	//计算从核接收数组数据长度和接收位置
	int len = cells / 64;
	int rest = cells % 64;
	int addr;
	if(CRTS_tid < rest){
		len++;
		addr = CRTS_tid * len;
	}else{
		addr = CRTS_tid * len + rest;
	}
	//接收数组数据
	CRTS_dma_iget(&p, slavePara.p + addr, len * sizeof(double), &DMARply);
	CRTS_dma_iget(&z, slavePara.z + addr, len * sizeof(double), &DMARply);
	DMARplyCount += 2;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
			
	//计算
	int i = 0;
	for(; i < len; i++){
		p[i] = z[i] + beta * p[i];
	}
	//传回计算结果
	CRTS_dma_iput(slavePara.p+addr, &p, len * sizeof(double), &DMARply);
	DMARplyCount++;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
}

#define SLICE_SIZE 1024

void csc_spmv_slave(CscChunkPara *chunk_para) {
    int cols; // chunk 的列数
    double *dst; // 计算结果应该放回主存的位置
    int slave_id = CRTS_smng_get_tid();

    // TODO: 根据 slave_id 对传入的 chunk_para 进行处理，获取到一些主存中的地址

    for (int i = 0; i < cols; ++i) {
        // 处理一列
        double result = .0;
        // TODO: 计算切片长度
        int n_slices;
        for (int j = 0; j < n_slices; ++j) {
            // 处理一个列的切片

        }
        CRTS_dma_iput(dst + i, &result, sizeof(double));
    }
    
    CRTS_sync_master_array();
}