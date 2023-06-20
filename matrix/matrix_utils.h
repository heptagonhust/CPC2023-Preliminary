#pragma once

#include <vector>

#include "pcg_def.h"

int count_ldu_matrix_nonzero_elements(const LduMatrix &ldu_matrix);
// void csr_to_csc(const CsrMatrix &csr_matrix, CscMatrix &csc_matrix);
void generate_ldu_matrix(LduMatrix &ldu_matrix, int size);
std::vector<double> plain_matrix_from_csr(const CsrMatrix &csr_matrix);
std::vector<double>
plain_matrix_from_splited_coo(const SplitedCooMatrix &splited);
// std::vector<double> plain_matrix_from_splited_csc_matrix(const SplitedCscMatrix &splited_matrix);
// std::vector<double> plain_matrix_from_csc(const CscMatrix &csc_matrix);
void free_ldu_matrix(LduMatrix &mtx);
