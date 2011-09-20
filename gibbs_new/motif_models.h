#ifndef MOTIF_MODELS_H
#define MOTIF_MODELS_H

#include "structs.h"

void initialize_motif_model(MotifModel *motif_model, int motif_width, int type);
int initialize_motif_models(int arg_position, char **argument);

void copy_motif_model(MotifModel* source, MotifModel* destination);

void free_motif_model(MotifModel* motif_model);
void free_motif_models();

void zero_counts(MotifModel *motif_model);

void increment_model_counts(MotifModel *motif_model, Sequence *sequence, Sample **sample, int motif_model_number);
void increment_all_counts(Sequence *sequence);

void decrement_counts(Sequence *sequence);

void update_motif_model_reverse_complement(MotifModel *forward_motif_model, MotifModel *reverse_motif_model);
void update_motif_model(MotifModel *motif_model);

void print_motif_model(MotifModel *motif_model);

#endif
