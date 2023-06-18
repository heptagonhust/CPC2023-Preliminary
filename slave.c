#include <slave.h>
#include "pcg_def.h"

#include <crts.h>
#include <swperf.h>

#include "ltypf.h"

typedef struct{
	double *p;
	double *z;
	double beta;
	int cells;
} Para;

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
//	penv_slave0_cycle_init();
//	penv_slave1_dcache_access_init();
//	penv_slave2_l1ic_access_init();
//	ltyperf_init();
	unsigned long icc1;	
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
//	penv_slave0_cycle_count(&icc1);
//	penv_slave1_dcache_access_count(&icc1);
//	penv_slave2_l1ic_access_count(&icc1);
//	ltyperf_report();
//	printf("%lu\n", icc1);



}
