#include <stdio.h>
#include <stdlib.h>

#include <algorithm>

#include "coo_matrix.h"
#include "pcg_def.h"

int main(void) {
    LduMatrix ldu;
    ldu.cells = 88362;
    ldu.faces = 263747;
    FILE *fp;

    ldu.diag = (double *)malloc(sizeof(double) * ldu.cells);
    ldu.lower = (double *)malloc(sizeof(double) * ldu.faces);
    ldu.upper = (double *)malloc(sizeof(double) * ldu.faces);
    ldu.lPtr = (int *)malloc(sizeof(int) * ldu.faces);
    ldu.uPtr = (int *)malloc(sizeof(int) * ldu.faces);

    fp = fopen("dia.txt", "r");
    for (int i = 0; i < ldu.cells; ++i) {
        fscanf(fp, "%lf", &ldu.diag[i]);
    }
    fp = fopen("lidx.txt", "r");
    for (int i = 0; i < ldu.faces; ++i) {
        fscanf(fp, "%d", &ldu.lPtr[i]);
    }
    fp = fopen("uidx.txt", "r");
    for (int i = 0; i < ldu.faces; ++i) {
        fscanf(fp, "%d", &ldu.uPtr[i]);
    }
    fp = fopen("lower.txt", "r");
    for (int i = 0; i < ldu.faces; ++i) {
        fscanf(fp, "%lf", &ldu.lower[i]);
    }
    fp = fopen("upper.txt", "r");
    for (int i = 0; i < ldu.faces; ++i) {
        fscanf(fp, "%lf", &ldu.upper[i]);
    }

    auto coo = ldu_to_splited_coo(ldu, 64);
    int max_col_size = 0, min_col_size = 999999999;
    int max_chunk_size = 0, min_chunk_size = 999999999;
    printf(
        "范围 col_num 非零元素数量 最大block大小 最小block大小 前四个block大小\n");
    for (int i = 0; i < coo.chunk_num; ++i) {
        auto &chunk = (*coo.chunks[i].chunk);
        int max_block_size = 0, min_block_size = 999999999;
        for (int j = 0; j < chunk.block_num; ++j) {
            int size =
                chunk.blocks[j + 1].block_off - chunk.blocks[j].block_off;
            max_block_size = std::max(max_block_size, size);
            min_block_size = std::min(min_block_size, size);
        }
        printf(
            "[%d,%d) %d %d %d %d",
            coo.chunk_ranges[i].col_begin,
            coo.chunk_ranges[i].col_begin + coo.chunk_ranges[i].col_num,
            coo.chunk_ranges[i].col_num,
            coo.chunk_ranges[i].size,
            max_block_size,
            min_block_size);
        for (int j = 0; j < 4; ++j) {
            int size =
                chunk.blocks[j + 1].block_off - chunk.blocks[j].block_off;
            printf(" %d", size);
        }
        putchar('\n');

        int size = coo.chunk_ranges[i].size;
        max_chunk_size = std::max(max_chunk_size, size);
        min_chunk_size = std::min(min_chunk_size, size);

        size = coo.chunk_ranges[i].col_num;
        max_col_size = std::max(max_col_size, size);
        min_col_size = std::min(min_col_size, size);
    }
    printf("最大chunk：%d\n", max_chunk_size);
    printf("最小chunk：%d\n", min_chunk_size);
    printf(
        "差距：%f\n",
        (max_chunk_size - min_chunk_size) * 100.0 / min_chunk_size);

    printf("最大col_num：%d\n", max_col_size);
    printf("最小col_num：%d\n", min_col_size);
    printf("差距：%f\n", (max_col_size - min_col_size) * 100.0 / min_col_size);
    return 0;
}
