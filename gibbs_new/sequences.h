#ifndef SEQUENCES_H
#define SEQUENCES_H

#include "structs.h"

void initialize_sequences_data();
void free_sequences();

long double foreground_probability(Sequence *sequence, MotifModel *motif_model, int end_position);

int possible_end_position(Sequence *sequence, MotifModel *motif_model, int position);

void calculate_site_counts(double *blocks);

void calculate_distance_to_invalid();

void calculate_background_site_probabilities(Sequence *sequence);
void calculate_site_probabilities(Sequence *sequence);

void read_sequence_file(char *filename);

void write_sites_to_file(FILE* file, const char *delimiter);


#endif
