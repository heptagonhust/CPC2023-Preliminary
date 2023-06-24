#include "spmv_master.h"
#include "swperf.h"
#include <crts.h>
#include <simd.h>
#define SLAVE_CORE_NUM 64


void spmv_para_from_splited_coo_matrix(
    const SplitedCooMatrix *mat,
    SpmvPara *para,
    int row_num,
    int col_num) {
    int chunk_num = mat->chunk_num;
    para->chunk_num = chunk_num;
    para->sp_row = row_num;
    para->sp_col = col_num;
    memcpy(para->chunks, mat->chunks, SLAVE_CORE_NUM * sizeof(SizedCooChunk));
    int max_row_num = 0;
    for (int i = 0; i < chunk_num; ++i) {
        int rows = para->chunks[i].chunk->row_end - para->chunks[i].chunk->row_begin;
        if (rows > max_row_num)
            max_row_num = rows;
    }
    para->max_block_row_num = max_row_num;
}
