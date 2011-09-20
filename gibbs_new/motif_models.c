#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "structs.h"
#include "motif_models.h"
#include "sequences.h"
#include "sampling.h"

extern int max_sites;

extern int number_motifs;
extern MotifModel **motif_models;

extern double foreground_pseudocounts;

void initialize_motif_model(MotifModel *motif_model, int motif_width, int type) {
    int j, k;

    motif_model->type = type;
    motif_model->motif_width = motif_width;

    motif_model->nucleotide_probabilities = (double**)malloc(sizeof(double*) * motif_width);
    
//    if (motif_model->type == MODEL_TYPE_REVERSE_COMPLEMENT) {
//        motif_model->reverse_nucleotide_probabilities = (double**)malloc(sizeof(double*) * motif_width);
//    }

    motif_model->counts = (double**)malloc(sizeof(double*) * motif_width);

    for (j = 0; j < motif_width; j++) {
        motif_model->counts[j] = (double*)malloc(sizeof(double) * ALPHABET_LENGTH);
        motif_model->nucleotide_probabilities[j] = (double*)malloc(sizeof(double) * ALPHABET_LENGTH);

//        if (motif_model->type == MODEL_TYPE_REVERSE_COMPLEMENT) {
//            motif_model->reverse_nucleotide_probabilities[j] = (double*)malloc(sizeof(double) * ALPHABET_LENGTH);
//        }

        for (k = 0; k < ALPHABET_LENGTH; k++) {
//            if (motif_model->type == MODEL_TYPE_REVERSE_COMPLEMENT) {
//                motif_model->reverse_nucleotide_probabilities[j][k] = 0.0;
//            }

            motif_model->nucleotide_probabilities[j][k] = 0.0;
            motif_model->counts[j][k] = 0.0;
        }
    }
}

int initialize_motif_models(int arg_position, char **arguments) {
    int i, j, k, motif_width;
    char type_string[100];
    char width_string[100];
    char tmp_char;
    int type;

    motif_models = (MotifModel**)malloc(sizeof(MotifModel*) * number_motifs);
    for (i = 0; i < number_motifs; i++) {
        motif_models[i] = (MotifModel*)malloc(sizeof(MotifModel));

        for (j = 0; j < (int)strlen(arguments[arg_position + i]); j++) {
            tmp_char = arguments[arg_position + i][j];

            if (tmp_char == ',') break;

            type_string[j] = tmp_char;
        }
        type_string[j] = '\0';
        fprintf(stderr, "type_string [%s]\n", type_string);

        for (k = j + 1; k < (int)strlen(arguments[arg_position + i]); k++) {
            tmp_char = arguments[arg_position + i][k];
            
            printf("tmp_char: %c\n", tmp_char);

            width_string[k - (j + 1)] = tmp_char;
        }
        width_string[k - (j + 1)] = '\0';
        fprintf(stderr, "width_string [%s]\n", width_string);
        
        motif_width = atoi(width_string);

        type = -1;

        if (!strcmp(type_string, "normal")) {
            type = MODEL_TYPE_NORMAL;

        } else if (!strcmp(type_string, "forward")) {
            type = MODEL_TYPE_FORWARD;

        } else if (!strcmp(type_string, "reverse")) {
            type = MODEL_TYPE_REVERSE;

        } else if (!strcmp(type_string, "palindromic")) {
            type = MODEL_TYPE_PALINDROMIC;

        } else {
            fprintf(stderr, "ERROR parsing motif model, unknown type: [%s] from argument [%s]\n", type_string, arguments[arg_position + i]);
            fprintf(stderr, "In file [%s], line [%d]\n", __FILE__, __LINE__);
            exit(0);
        }

        initialize_motif_model(motif_models[i], motif_width, type);
    }

    return number_motifs - 1;
}

void copy_motif_model(MotifModel* source, MotifModel* destination) {
    int j, k;

    if (source->motif_width != destination->motif_width) {
        fprintf(stderr, "ERROR copying motif model, source->width [%d] != destination->width [%d]\n", source->motif_width, destination->motif_width);
        fprintf(stderr, "In file [%s], line [%d]\n", __FILE__, __LINE__);
        exit(0);
    }

    for (j = 0; j < source->motif_width; j++) {
        for (k = 0; k < ALPHABET_LENGTH; k++) {
            destination->nucleotide_probabilities[j][k] = source->nucleotide_probabilities[j][k];
            destination->counts[j][k] = source->counts[j][k];
        }
    }
}

void free_motif_model(MotifModel* motif_model) {
    int j;
    for (j = 0; j < motif_model->motif_width; j++) {
        free(motif_model->counts[j]);
        free(motif_model->nucleotide_probabilities[j]);
    }
}

void free_motif_models() {
    int i;

    for (i = 0; i < number_motifs; i++) {
        free_motif_model(motif_models[i]);
        free(motif_models[i]);
    }
    free(motif_models);
}

void zero_counts(MotifModel *motif_model) {
    int j;

    for (j = 0; j < motif_model->motif_width; j++) {
        motif_model->counts[j][A_POSITION] = 0;
        motif_model->counts[j][C_POSITION] = 0;
        motif_model->counts[j][G_POSITION] = 0;
        motif_model->counts[j][T_POSITION] = 0;
    }
}

void increment_counts(Sequence *sequence, MotifModel *motif_model, Sample *sample) {
    int position, end_position, j;
    char nucleotide;

    end_position = sample->end_position + 1;

    if (end_position - motif_model->motif_width < 0) {
        fprintf(stderr, "ERROR (increment_counts): trying to read character with position < 0\n");

        fprintf(stderr, "subsequence examined [");
        for (j = 0; j < motif_model->motif_width; j++) {
            fprintf(stderr, "%c", sequence->nucleotides[end_position - motif_model->motif_width + j]);
        }
        fprintf(stderr, "]\n");

        exit(0);
    }

    for (position = 0; position < motif_model->motif_width; position++) {
        nucleotide = sequence->nucleotides[end_position - motif_model->motif_width + position];

        switch (nucleotide) {
            case 'A': motif_model->counts[position][A_POSITION]++; break;
            case 'C': motif_model->counts[position][C_POSITION]++; break;
            case 'G': motif_model->counts[position][G_POSITION]++; break;
            case 'T': motif_model->counts[position][T_POSITION]++; break;

            default:
                      fprintf(stderr, "ERROR (increment_model_counts): unknown character found in sequence [%s][%s][%d] at position [%d] nucleotide position [%d] -- [%c] [%d]\n", sequence->name, sequence->nucleotides, sequence->length, position, (end_position - motif_model->motif_width + position), nucleotide, (int)nucleotide);
                      fprintf(stderr, "sample->end_position [%d], motif_width [%d]\n", sample->end_position, motif_model->motif_width);

                      fprintf(stderr, "subsequence examined [");
                      for (j = 0; j < motif_model->motif_width; j++) {
                          fprintf(stderr, "%c", sequence->nucleotides[end_position - motif_model->motif_width + j]);
                      }
                      fprintf(stderr, "]\n");
                      fprintf(stderr, "Error in file [%s], line [%d]\n", __FILE__, __LINE__);
        }
    }
}

void increment_model_counts(MotifModel *motif_model, Sequence *sequence, Sample **sample, int motif_model_number) {
    int i;
    
    for (i = 0; i < max_sites; i++) {
        if (sample[i]->motif_model != motif_model_number) continue;

        increment_counts(sequence, motif_model, sample[i]);
    }
}

void increment_all_counts(Sequence *sequence) {
    int i;
    MotifModel *motif_model;
    Sample *sample;

    for (i = 0; i < max_sites; i++) {
        sample = sequence->sampled_sites[i];

        if (sample->motif_model < 0 || sample->end_position < 0) continue;

        motif_model = motif_models[ sample->motif_model ];

        increment_counts(sequence, motif_model, sample);
    }
}

void decrement_counts(Sequence *sequence) {
    char nucleotide;
    int i, j;
    int position;
    int end_position;
    MotifModel *motif_model;
    Sample *sample;

    for (i = 0; i < max_sites; i++) {
        sample = sequence->sampled_sites[i];

        if (sample->motif_model < 0 || sample->end_position < 0) continue;
        motif_model = motif_models[ sample->motif_model ];
        end_position = sample->end_position + 1;

        if (end_position - motif_model->motif_width < 0) {
            fprintf(stderr, "ERROR (decrement_counts): trying to read character with position < 0\n");

            fprintf(stderr, "subsequence examined [");
            for (j = 0; j < motif_model->motif_width; j++) {
                fprintf(stderr, "%c", sequence->nucleotides[end_position - motif_model->motif_width + j]);
            }
            fprintf(stderr, "]\n");

            exit(0);
        }

        for (position = 0; position < motif_model->motif_width; position++) {
            nucleotide = sequence->nucleotides[end_position - motif_model->motif_width + position];

    //        printf("incrementing count for position [%d] for nucleotide [%c]\n", position, nucleotide);

            switch (nucleotide) {
                case 'A': motif_model->counts[position][A_POSITION]--; break;
                case 'C': motif_model->counts[position][C_POSITION]--; break;
                case 'G': motif_model->counts[position][G_POSITION]--; break;
                case 'T': motif_model->counts[position][T_POSITION]--; break;

                default:
                          fprintf(stderr, "ERROR (decrement_counts): unknown character found in sequence [%s][%s][%d] at position [%d] nucleotide position [%d] -- [%c] [%d]\n", sequence->name, sequence->nucleotides, sequence->length, position, (end_position - motif_model->motif_width + position), nucleotide, (int)nucleotide);
                          fprintf(stderr, "sample->end_position [%d], motif_width [%d]\n", sample->end_position, motif_model->motif_width);

                          fprintf(stderr, "this a possible end position? %d\n", possible_end_position(sequence, motif_model, position));

                          fprintf(stderr, "subsequence examined [");
                          for (j = 0; j < motif_model->motif_width; j++) {
                              fprintf(stderr, "%c", sequence->nucleotides[end_position - motif_model->motif_width + j]);
                          }
                          fprintf(stderr, "]\n");

            }
        }
    }
}

void update_motif_model_reverse_complement(MotifModel *forward_motif_model, MotifModel *reverse_motif_model) {
    double count_a, count_c, count_g, count_t;
    double sum;
    int position;

    switch (forward_motif_model->type) {
        case MODEL_TYPE_FORWARD:
                for (position = 0; position < forward_motif_model->motif_width; position++) {

                    count_a = foreground_pseudocounts + forward_motif_model->counts[position][A_POSITION] + reverse_motif_model->counts[forward_motif_model->motif_width - (position + 1)][A_POSITION_REVERSE];
                    count_c = foreground_pseudocounts + forward_motif_model->counts[position][C_POSITION] + reverse_motif_model->counts[forward_motif_model->motif_width - (position + 1)][C_POSITION_REVERSE];
                    count_g = foreground_pseudocounts + forward_motif_model->counts[position][G_POSITION] + reverse_motif_model->counts[forward_motif_model->motif_width - (position + 1)][G_POSITION_REVERSE];
                    count_t = foreground_pseudocounts + forward_motif_model->counts[position][T_POSITION] + reverse_motif_model->counts[forward_motif_model->motif_width - (position + 1)][T_POSITION_REVERSE];

                    sum = count_a + count_c + count_g + count_t;

                    forward_motif_model->nucleotide_probabilities[position][A_POSITION] = count_a / sum; 
                    forward_motif_model->nucleotide_probabilities[position][C_POSITION] = count_c / sum; 
                    forward_motif_model->nucleotide_probabilities[position][G_POSITION] = count_g / sum; 
                    forward_motif_model->nucleotide_probabilities[position][T_POSITION] = count_t / sum; 

                    reverse_motif_model->nucleotide_probabilities[forward_motif_model->motif_width - (position + 1)][A_POSITION_REVERSE] = forward_motif_model->nucleotide_probabilities[position][A_POSITION];
                    reverse_motif_model->nucleotide_probabilities[forward_motif_model->motif_width - (position + 1)][C_POSITION_REVERSE] = forward_motif_model->nucleotide_probabilities[position][C_POSITION];
                    reverse_motif_model->nucleotide_probabilities[forward_motif_model->motif_width - (position + 1)][G_POSITION_REVERSE] = forward_motif_model->nucleotide_probabilities[position][G_POSITION];
                    reverse_motif_model->nucleotide_probabilities[forward_motif_model->motif_width - (position + 1)][T_POSITION_REVERSE] = forward_motif_model->nucleotide_probabilities[position][T_POSITION];
                }
                break;

        default:
                fprintf(stderr, "ERROR: unknown model type in update_motif_model_reverse_complement, FILE [%s], LINE [%d]\n", __FILE__, __LINE__);
                fprintf(stderr, "FORWARD MODEL:\n");
                print_motif_model(forward_motif_model);

                fprintf(stderr, "\n");
                fprintf(stderr, "REVERSEMODEL:\n");
                print_motif_model(reverse_motif_model);
                exit(0);
    }

}

void update_motif_model(MotifModel *motif_model) {
    double count_a, count_c, count_g, count_t;
    double sum;
    int position;

    switch (motif_model->type) {
        case MODEL_TYPE_NORMAL:
                for (position = 0; position < motif_model->motif_width; position++) {

                    count_a = foreground_pseudocounts + motif_model->counts[position][A_POSITION];
                    count_c = foreground_pseudocounts + motif_model->counts[position][C_POSITION];
                    count_g = foreground_pseudocounts + motif_model->counts[position][G_POSITION];
                    count_t = foreground_pseudocounts + motif_model->counts[position][T_POSITION];

                    sum = count_a + count_c + count_g + count_t;

                    motif_model->nucleotide_probabilities[position][A_POSITION] = count_a / sum; 
                    motif_model->nucleotide_probabilities[position][C_POSITION] = count_c / sum; 
                    motif_model->nucleotide_probabilities[position][G_POSITION] = count_g / sum; 
                    motif_model->nucleotide_probabilities[position][T_POSITION] = count_t / sum; 

            //        printf("set probs[%d] to %10.5lf %10.5lf %10.5lf %10.5lf\n", position, (count_a + 0.28) / sum, (count_c + 0.28) / sum, (count_g + 0.28) / sum, (count_t + 0.28) / sum);
                }
                break;

        case MODEL_TYPE_PALINDROMIC:
                for (position = 0; position < (motif_model->motif_width / 2); position++) {
//                    printf("calculating palindromic probability for position [%d] and [%d]\n", position, (motif_model->motif_width - (position + 1)));

//                    printf("foreground_pseudocounts: %lf\n", foreground_pseudocounts);
                    count_a = foreground_pseudocounts + motif_model->counts[position][A_POSITION] + motif_model->counts[motif_model->motif_width - (position + 1)][A_POSITION_REVERSE];
                    count_c = foreground_pseudocounts + motif_model->counts[position][C_POSITION] + motif_model->counts[motif_model->motif_width - (position + 1)][C_POSITION_REVERSE];
                    count_g = foreground_pseudocounts + motif_model->counts[position][G_POSITION] + motif_model->counts[motif_model->motif_width - (position + 1)][G_POSITION_REVERSE];
                    count_t = foreground_pseudocounts + motif_model->counts[position][T_POSITION] + motif_model->counts[motif_model->motif_width - (position + 1)][T_POSITION_REVERSE];

                    sum = count_a + count_c + count_g + count_t;

                    motif_model->nucleotide_probabilities[position][A_POSITION] = count_a / sum; 
                    motif_model->nucleotide_probabilities[position][C_POSITION] = count_c / sum; 
                    motif_model->nucleotide_probabilities[position][G_POSITION] = count_g / sum; 
                    motif_model->nucleotide_probabilities[position][T_POSITION] = count_t / sum; 

                    motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][T_POSITION_REVERSE] = count_t / sum; 
                    motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][G_POSITION_REVERSE] = count_g / sum; 
                    motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][C_POSITION_REVERSE] = count_c / sum; 
                    motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][A_POSITION_REVERSE] = count_a / sum; 

////                    printf("probabilities: [%lf][%lf][%lf][%lf]\n", motif_model->nucleotide_probabilities[position][A_POSITION], motif_model->nucleotide_probabilities[position][C_POSITION], motif_model->nucleotide_probabilities[position][G_POSITION], motif_model->nucleotide_probabilities[position][T_POSITION]);

//                    printf("probabilities: [%lf][%lf][%lf][%lf]\n", motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][A_POSITION], motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][C_POSITION], motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][G_POSITION], motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][T_POSITION]);
                }

                if ((motif_model->motif_width % 2) == 1) {
                    position = (motif_model->motif_width / 2) + 1;
//                    printf("palindromic middle for model width [%d] is [%d]\n", motif_model->motif_width, position);

                    count_a = foreground_pseudocounts + motif_model->counts[position][A_POSITION];
                    count_c = foreground_pseudocounts + motif_model->counts[position][C_POSITION];
                    count_g = foreground_pseudocounts + motif_model->counts[position][G_POSITION];
                    count_t = foreground_pseudocounts + motif_model->counts[position][T_POSITION];

                    sum = count_a + count_c + count_g + count_t;

                    motif_model->nucleotide_probabilities[position][A_POSITION] = count_a / sum; 
                    motif_model->nucleotide_probabilities[position][C_POSITION] = count_c / sum; 
                    motif_model->nucleotide_probabilities[position][G_POSITION] = count_g / sum; 
                    motif_model->nucleotide_probabilities[position][T_POSITION] = count_t / sum; 
                }

                break;
        default:
                fprintf(stderr, "ERROR: unknown model type in update_motif_model, FILE [%s], LINE [%d]\n", __FILE__, __LINE__);
                print_motif_model(motif_model);
                exit(0);
    }
}

void print_motif_model(MotifModel *motif_model) {
    int j, k;

    printf("Model type: %d\n", motif_model->type);
    printf("Model length: %d\n", motif_model->motif_width);

    for (j = 0; j < ALPHABET_LENGTH; j++) {
        printf("letter count [%c]: ", ALPHABET[j]);

        for (k = 0; k < motif_model->motif_width; k++) {
            printf("%10.5lf    ", motif_model->counts[k][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (j = 0; j < ALPHABET_LENGTH; j++) {
        printf("letter prob [%c]: ", ALPHABET[j]);

        for (k = 0; k < motif_model->motif_width; k++) {
            printf("%5.10lf    ", motif_model->nucleotide_probabilities[k][j]);
        }
        printf("\n");
    }
    
    printf("\n");

}




