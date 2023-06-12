#include "spmv.h"

#include <crts.h>
#include <slave.h>
#include "slave_def.h"

#include "spmv_def.h"

#include "slave_def.h"

__thread_local crts_rply_t DMARply = 0;
__thread_local unsigned int DMARplyCount = 0;

#define BUF_ITEM 1024
__thread_local double redu_buf[BUF_ITEM] __attribute__((aligned(64)));

void csc_spmv_slave(void *para) {
    SpmvPara spmv_para;
    DMA_GET(&spmv_para, para, sizeof(SpmvPara), &DMARply);
    int id = CRTS_tid;

    int col = spmv_para.sp_col;
    int row = spmv_para.sp_row;

    double *chunk_result = CRTS_pldm_malloc(sizeof(double) * row);
    for (double *p = chunk_result; p < chunk_result + row; ++p) {
        *p = 0;
    }

    // 获取 chunk
    CscChunk **chunk_array_ptr = spmv_para.chunks;
    int chunk_num = spmv_para.chunk_num;
    CscChunk **chunk_array = CRTS_pldm_malloc(sizeof(CscChunk *) * chunk_num);
    DMA_GET(chunk_array, chunk_array_ptr, sizeof(CscChunk *) * chunk_num, &DMARply);
    CscChunk *chunk_ptr = chunk_array[id];
    CRTS_pldm_free(chunk_array, sizeof(CscChunk *) * chunk_num);

    int size;
    DMA_GET(&size, chunk_ptr, sizeof(int), &DMARply);
    CscChunk *chunk = CRTS_pldm_malloc(size);
    DMA_GET(chunk, chunk_ptr, size, &DMARply);
    int block_num = chunk->block_num;

    // 获取该 chunk 对应 vec 的 slice
    double *vec = chunk->vec;
    int col_begin = chunk->col_begin;
    int col_end = chunk->col_end;
    int slice_size = col_end - col_begin;
    double *slice = CRTS_pldm_malloc(slice_size);
    DMA_GET(slice, vec + col_begin, slice_size, &DMARply);

    double *block_result = chunk_result;
    CscBlock block;
    int buff_size = 0;
    void *block_buff = 0;
    for (int i = 0; i < block_num; ++i) {
        CscBlock *block_ptr = &chunk->blocks[i];
        DMA_GET(&block, block_ptr, sizeof(CscBlock), &DMARply);
        int col_num = block.col_num;
        int data_size = block.data_size;
        int block_size = (sizeof(double) + sizeof(int)) * data_size
            + sizeof(int) * (block.col_num + 1);

        // 每次比较一下，减少分配次数
        if (block_size > buff_size) {
            if (block_buff != 0) {
                CRTS_pldm_free(block_buff, block_size);
            }
            block_buff = CRTS_pldm_malloc(block_size);
            buff_size = block_size;
        }

        // TODO: 访存优化，需要优化这三个数组的排布顺序
        double *data = (double *)block_buff;
        int *rows = (int *)(data + data_size);
        int *col_off = rows + data_size;

        double *data_ptr = block.data;
        int *rows_ptr = block.rows;
        int *col_off_ptr = block.col_off;

        DMA_IGET(data, data_ptr, sizeof(double) * data_size, &DMARply);
        DMA_IGET(rows, rows_ptr, sizeof(double) * data_size, &DMARply);
        DMA_IGET(col_off, col_off_ptr, sizeof(double) * data_size, &DMARply);
        DMA_WAIT;

        // block 与 slice 的 spmv
        for (int j = 0; j < col_num; ++j) {
            int start = block.col_off[j];
            int end = block.col_off[j + 1];
            int num = end - start;
            double tmp = 0.;
            for (int k = 0; k < num; ++k) {
                tmp += slice[rows[start + k]] * data[start + k];
            }
            block_result[j] = tmp;
        }

        block_result += col_num;
    }

    CRTS_pldm_free(block_buff, buff_size);

    // 规约操作
    CRTS_scoll_redurt(
        chunk_result,
        chunk_result,
        row,
        CRTS_double,
        OP_add,
        redu_buf,
        BUF_ITEM);

    // slave core 0 将结果写回主存
    if (id == 0) {
        double *result = spmv_para.result;
        DMA_PUT(result, chunk_result, sizeof(double) * row, &DMARply);
    }

    CRTS_pldm_free(chunk_result, sizeof(double) * row);
    CRTS_pldm_free(slice, slice_size);
    CRTS_pldm_free(chunk, size);
}