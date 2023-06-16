#pragma once

#include "pcg.h"
#include "spmv_def.h"

extern void ldu_to_csc(const LduMatrix &ldu_matrix, CscMatrix &csc_matrix);
extern void free_csc_matrix(CscMatrix &mtx);
extern void csc_chunk_pack(CscChunk &chunk);
extern void csc_chunk_unpack(CscChunk &chunk);
extern void free_packed_splited_csc_matrix_chunk(CscChunk *chunk);
extern void free_packed_splited_csc_matrix(SplitedCscMatrix *splited);
extern SplitedCscMatrix split_csc_matrix(const CscMatrix &mtx, int chunk_num);
