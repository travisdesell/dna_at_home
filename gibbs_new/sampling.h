#ifndef SAMPLING_H
#define SAMPLING_H

#include "structs.h"

void sample_uniform_random(Sequence *sequence);

void resample_from_models(Sequence *sequence);

void print_sample_and_nearest(Sequence *sequence, MotifModel *motif_model, int end_position, int nearest);

void print_sample(Sequence *sequence, MotifModel *motif_model, int end_position);

#endif
