#ifndef SHIFTING_H
#define SHIFTING_H

#include "structs.hpp"

void initialize_shifting(vector<MotifModel> &shifted_motif_models);

void attempt_shifting(int max_shift_distance, vector<Sequence> &sequences, vector<MotifModel> &motif_models);

#endif
