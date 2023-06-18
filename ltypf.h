// #ifndef _LTYPF_H_
// #define _LTYPF_H_

#include <stdio.h>
#include <swperf.h>
// init
void ltyperf_init(void) {
    penv_slave0_cycle_init();// cycle passed
    penv_slave1_dcache_access_init();// dcache access count
    penv_slave2_l1ic_access_init();// l1 instruction cache 
    penv_slave3_dcache_hit_ldm_init(); // dcache hit ldm
    penv_slave4_readldm_dma_rma_init(); //ldm dma rma count
    penv_slave5_dcache_access_ldm_init(); 

}

// report
void ltyperf_report(void) {
    unsigned long _answer0 = 0;
    penv_slave0_cycle_count(&_answer0);
    printf("slave cycle passed: %lu cycles \n", _answer0);

    unsigned long _answer1 = 0;
    penv_slave1_dcache_access_count(&_answer1);
    printf("dcache access count: %lu counts \n", _answer1);
 
    unsigned long _answer2 = 0;
    penv_slave2_l1ic_access_count(&_answer2);
    printf("l1 instruction cache access count: %lu counts \n", _answer2);

    unsigned long _answer3 = 0;	
    penv_slave3_dcache_hit_ldm_init(&_answer3);
    printf("dcache ldm hit count: %lu counts \n", _answer3);

    unsigned long _answer4 = 0;
    penv_slave4_readldm_dma_rma_count(&_answer4);
    printf("ldm rma dma access count: %lu counts \n", _answer4);
   
    unsigned long _answer5 = 0;	
    penv_slave5_dcache_access_ldm_count(&_answer5);
    printf("dcache access ldm count: %lu counts \n", _answer5);
 
    unsigned long _answer6 = 0;                      
    unsigned long _answer7 = 0;

}        	
