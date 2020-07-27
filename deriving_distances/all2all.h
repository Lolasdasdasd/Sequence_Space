#ifndef ALL2ALL_H
#define ALL2ALL_H

#include <cstddef>
#include <string>
#include "model.h"


pair<uint64_t, uint64_t> fill_one_block(uint64_t pos, uint64_t pos2, size_t lf, size_t domain_size, Model* model, StringSet<Peptide>& block1, StringSet<Peptide>& block2, uint64_t block_index, uint64_t block_stop);

void all2all(size_t lf, size_t threads, size_t domain_size,  size_t start_seed, std::string outfile, std::string path, Model* model);

#endif



