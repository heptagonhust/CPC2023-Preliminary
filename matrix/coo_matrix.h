#pragma once

#include "pcg.h"
#include "spmv_def.h"

extern SplitedCooMatrix
ldu_to_splited_coo(const LduMatrix &ldu_matrix, int chunk_num);
extern void
free_splited_coo(SplitedCooMatrix &);
