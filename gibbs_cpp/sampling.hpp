#ifndef SAMPLING_H
#define SAMPLING_H

#include "sequences.hpp"

using namespace std;

void print_sample_and_nearest(Sequence *sequence, MotifModel &motif_model, int end_position, int nearest);

void print_sample(Sequence *sequence, MotifModel &motif_model, int end_position);

#endif
