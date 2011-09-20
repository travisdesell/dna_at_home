#ifdef _BOINC_
    #include "diagnostics.h"
    #include "util.h"
    #include "filesys.h"
    #include "boinc_api.h"
    #include "mfile.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "structs.h"
#include "sequences.h"
#include "motif_models.h"

extern int max_sites;

extern int number_sequences;
extern Sequence **sequences;

extern int number_motifs;
extern MotifModel **motif_models;

extern int max_shift_distance;

extern double background_pseudocounts;
extern double background_nucleotide_probability[ALPHABET_LENGTH];

void initialize_sequences_data() {
    int i, j, k;
    Sequence *current_sequence;

    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];

        current_sequence->possible_sites = (long double*)malloc(sizeof(long double) * (max_sites + 1));
        current_sequence->possible_site_counts = (long double**)malloc(sizeof(long double*) * max_sites);

        current_sequence->background_site_probability = (long double**)malloc(sizeof(long double*) * number_motifs);
        current_sequence->site_probability_ratio = (long double**)malloc(sizeof(long double*) * number_motifs);

        current_sequence->site_probability = (long double**)malloc(sizeof(long double*) * max_sites);
        current_sequence->motif_probability_contribution = (long double***)malloc(sizeof(long double**) * max_sites);

        current_sequence->sampled_sites = (Sample**)malloc(sizeof(Sample*) * max_sites);

        for (j = 0; j < max_sites; j++) {
            current_sequence->sampled_sites[j] = (Sample*)malloc(sizeof(Sample));

            current_sequence->possible_site_counts[j] = (long double*)malloc(sizeof(long double) * current_sequence->length);
            current_sequence->site_probability[j] = (long double*)malloc(sizeof(long double) * current_sequence->length);
            current_sequence->motif_probability_contribution[j] = (long double**)malloc(sizeof(long double*) * current_sequence->length);

            for (k = 0; k < current_sequence->length; k++) {
                current_sequence->motif_probability_contribution[j][k] = (long double*)malloc(sizeof(long double) * number_motifs);
            }
        }

        for (j = 0; j < number_motifs; j++) {
            current_sequence->background_site_probability[j] = (long double*)malloc(sizeof(long double) * current_sequence->length);
            current_sequence->site_probability_ratio[j] = (long double*)malloc(sizeof(long double) * current_sequence->length);
        }

        /**
         *  Only needed if we are doing shifting.
         */
//        printf("max_shift_distance: %d\n", max_shift_distance);
        if (max_shift_distance > 0) {
//            printf("mallocing shifted sites\n");
            current_sequence->shifted_sites = (Sample**)malloc(sizeof(Sample*) * max_sites);

            for (j = 0; j < max_sites; j++) {
                current_sequence->shifted_sites[j] = (Sample*)malloc(sizeof(Sample));
            }
        }
    }
}

void free_sequences() {
    int i, j, k;
    Sequence *current_sequence;

    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];
        free(current_sequence->possible_sites);

        for (j = 0; j < max_sites; j++) {
            for (k = 0; k < current_sequence->length; k++) {
                free(current_sequence->motif_probability_contribution[j][k]);
            }
            free(current_sequence->possible_site_counts[j]);
            free(current_sequence->site_probability[j]);
            free(current_sequence->motif_probability_contribution[j]);

            free(current_sequence->sampled_sites[j]);
        }
        free(current_sequence->possible_site_counts);
        free(current_sequence->site_probability);
        free(current_sequence->motif_probability_contribution);

        for (j = 0; j < number_motifs; j++) {
            free(current_sequence->background_site_probability[j]);
            free(current_sequence->site_probability_ratio[j]);
        }
        free(current_sequence->background_site_probability);
        free(current_sequence->site_probability_ratio);

        if (max_shift_distance > 0) {
            for (j = 0; j < max_sites; j++) {
                free(current_sequence->shifted_sites[j]);
            }
            free(current_sequence->shifted_sites);
        }
    }

    free(sequences);
}


int is_valid_nucleotide(char nucleotide) {
    switch (nucleotide) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return 1;
        default:
            return 0;
    }
}

int get_nucleotide_position(char nucleotide) {
    switch (nucleotide) {
        case 'A': return A_POSITION;
        case 'C': return C_POSITION;
        case 'G': return G_POSITION;
        case 'T': return T_POSITION;
        default: return -1;
    }
}

void sequences_to_uppercase() {
    int i, j;
    Sequence *current_sequence;

    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];

        for (j = 0; j < current_sequence->length; j++) {
            current_sequence->nucleotides[j] = toupper(current_sequence->nucleotides[j]);
        }
    }
}

void calculate_background_nucleotide_probabilities() {
    int i, j;
    int total, a_count, c_count, g_count, t_count;
    Sequence *current_sequence;

    total = background_pseudocounts + background_pseudocounts +  background_pseudocounts + background_pseudocounts;
    a_count = background_pseudocounts;
    c_count = background_pseudocounts;
    g_count = background_pseudocounts;
    t_count = background_pseudocounts;

    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];

        for (j = 0; j < current_sequence->length; j++) {
            switch (current_sequence->nucleotides[j]) {
                case 'A': a_count++;
                          total++;
                          break;

                case 'C': c_count++;
                          total++;
                          break;

                case 'G': g_count++;
                          total++;
                          break;
                          
                case 'T': t_count++;
                          total++;
                          break;

//                default:
//                          fprintf(stderr, "Unknown letter in sequence [%d] at position [%d]: '%c'\n", i, j, current_sequence->nucleotides[j]);
//                          exit(0);
            }
        }


    }

    background_nucleotide_probability[A_POSITION] = (double)(a_count + t_count) / (double)total;
    background_nucleotide_probability[C_POSITION] = (double)(c_count + g_count) / (double)total;
    background_nucleotide_probability[G_POSITION] = (double)(c_count + g_count) / (double)total;
    background_nucleotide_probability[T_POSITION] = (double)(a_count + t_count) / (double)total;
//    background_nucleotide_probability[A_POSITION] = (double)(a_count) / (double)total;
//    background_nucleotide_probability[C_POSITION] = (double)(c_count) / (double)total;
//    background_nucleotide_probability[G_POSITION] = (double)(g_count) / (double)total;
//    background_nucleotide_probability[T_POSITION] = (double)(t_count) / (double)total;

    printf("background probabilities:\n");
    printf("A: %lf\n", background_nucleotide_probability[A_POSITION]);
    printf("C: %lf\n", background_nucleotide_probability[C_POSITION]);
    printf("G: %lf\n", background_nucleotide_probability[G_POSITION]);
    printf("T: %lf\n", background_nucleotide_probability[T_POSITION]);
}

void calculate_distance_to_invalid() {
    int i, j;
    Sequence *current_sequence;

    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];
        current_sequence->distance_to_previous_invalid = (int*)malloc(sizeof(int) * current_sequence->length);

        if (is_valid_nucleotide(current_sequence->nucleotides[0])) {
            current_sequence->distance_to_previous_invalid[0] = 1;
        } else {
            current_sequence->distance_to_previous_invalid[0] = 0;
        }

        for (j = 1; j < current_sequence->length; j++) {
            if (is_valid_nucleotide(current_sequence->nucleotides[j])) {
                current_sequence->distance_to_previous_invalid[j] = current_sequence->distance_to_previous_invalid[j-1] + 1;
            } else {
                current_sequence->distance_to_previous_invalid[j] = 0;
            }
        }
    }
}

int possible_end_position(Sequence *sequence, MotifModel *motif_model, int position) {
    if (position < 1 || position >= sequence->length || sequence->distance_to_previous_invalid[position] < motif_model->motif_width) return 0;
    return 1;
}


/**
 *  This counts the number of possible ways to generate 0...max_sites within a sequence
 */
void calculate_site_counts(double *blocks) {
    int i, j, k, l;
    Sequence *current_sequence;
    MotifModel *current_motif_model;

//    printf("calculating counts\n");

    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];

        for (j = 0; j < max_sites; j++) {

            for (k = 0; k < current_sequence->length; k++) {
//                printf("number_motifs [%d], current_sequence->length [%d]\n", number_motifs, current_sequence->length);

                if (k == 0) {
                    current_sequence->possible_site_counts[j][0] = 0;
                } else {
                    current_sequence->possible_site_counts[j][k] = current_sequence->possible_site_counts[j][k - 1];
                }

                for (l = 0; l < number_motifs; l++) {
                    current_motif_model = motif_models[l];

                    if (possible_end_position(current_sequence, current_motif_model, k)) {
                        //j == 0 is a special case, as all counts[-1][k] are == 1 (where -1 is 0 sites)
                        if (j == 0) {
                            if (k - current_motif_model->motif_width + 1 >= 0) {
                                current_sequence->possible_site_counts[j][k] += 1;
                            }
                        } else if (k - current_motif_model->motif_width >= 0) {
                            current_sequence->possible_site_counts[j][k] += current_sequence->possible_site_counts[j - 1][k - current_motif_model->motif_width];
                        }
                    }
                    
//                    printf("possible_site_counts[%d][%d]: %Lf\n", j, k, current_sequence->possible_site_counts[j][k]);
                }
            }

        }

        /*
        printf("sequence: [%s]\n", current_sequence->nucleotides);
        for (j = 0; j < max_sites; j++) {
            printf("number_sites [%d]:", j);
            for (k = 0; k < current_sequence->length; k++) {
                printf(" %10.2Lf", current_sequence->possible_site_counts[j][k]);
            }
            printf("\n");
        }
        */


        current_sequence->possible_sites[0] = 1 * blocks[0]; //theres one way to have 0 sites sampled -- but we multiply this by blocks
        current_sequence->total_possible_sites = 0;

//        printf("sequence [%s]\n", current_sequence->nucleotides);
//        printf("possible_sites[0]: 1.0\n");

        for (j = 0; j < max_sites; j++) {
            current_sequence->possible_sites[j+1] = blocks[j+1] * current_sequence->possible_site_counts[j][current_sequence->length - 1];
        }

        for (j = 0; j < max_sites + 1; j++) {
            current_sequence->total_possible_sites += current_sequence->possible_sites[j];
//            printf("possible_sites[%d]: %Lf\n", j, current_sequence->possible_sites[j]);
        }
//        printf("total_possible_sites: %Lf\n", current_sequence->total_possible_sites);
    }
}

long double background_probability(Sequence *sequence, double *background_nucleotide_probability, MotifModel *motif_model, int end_position) {
    int i;
    int position;
    int nucleotide_position;
    long double probability;

    probability = 1.0;

    end_position++;

    for (i = 0; i < motif_model->motif_width; i++) {
        position = end_position - motif_model->motif_width + i;
        nucleotide_position = get_nucleotide_position( sequence->nucleotides[position] );

        probability *= background_nucleotide_probability[nucleotide_position];
    }

    return probability;
}

long double foreground_probability(Sequence *sequence, MotifModel *motif_model, int end_position) {
    int i;
    int position;
    int nucleotide_position;
    long double probability, reverse_probability;

    probability = 1.0;
    reverse_probability = 1.0;

    end_position++;
    for (i = 0; i < motif_model->motif_width; i++) {
        position = end_position - motif_model->motif_width + i;
        nucleotide_position = get_nucleotide_position( sequence->nucleotides[position] );

        probability *= motif_model->nucleotide_probabilities[i][nucleotide_position];
    }

    if (motif_model->type == MODEL_TYPE_FORWARD || motif_model->type == MODEL_TYPE_REVERSE) probability /= 2.0;

    return probability;
}

void calculate_background_site_probabilities(Sequence *sequence) {
    int i, j;
    MotifModel *current_motif;

    /**
     * Currently there's no reason we need to recalculate the background_probabilities -- this would only be done if the motifs are resized
     */
    for (i = 0; i < number_motifs; i++) {
        current_motif = motif_models[i];

//        print_motif_model(current_motif);

        for (j = 0; j < sequence->length; j++) {
            if (possible_end_position(sequence, current_motif, j)) {
                sequence->background_site_probability[i][j] = background_probability(sequence, background_nucleotide_probability, current_motif, j);
            } else {
                sequence->background_site_probability[i][j] = 0;
            }
//          printf("%.20Lf ", sequence->background_site_probability[i][j]);
        }
    }
}

void calculate_site_probabilities(Sequence *sequence) {
    int i, j, k, l;
    long double motif_contribution;
    MotifModel *current_motif;

//    printf("calculated background probabilities\n");

    /**
     *  We always need to recalculate the site probability ratios (forground / background)
     */
//    printf("caculating site probability ratios: \n");
    for (i = 0; i < number_motifs; i++) {
        current_motif = motif_models[i];

        for (j = 0; j < sequence->length; j++) {
            if (possible_end_position(sequence, current_motif, j)) {
                sequence->site_probability_ratio[i][j] = foreground_probability(sequence, current_motif, j) / sequence->background_site_probability[i][j];

                if (sequence->site_probability_ratio[i][j] < 0) {
                    fprintf(stderr, "ERROR: calculated site probability ratio < 0: file [%s] line [%d]\n", __FILE__, __LINE__);
                    fprintf(stderr, "sequence->site_probability_ratio[%d][%d]: %.30Lf\n", i, j, sequence->site_probability_ratio[i][j]);
                    fprintf(stderr, "foreground_prob: %.30Lf, background_prob: %.30Lf\n", foreground_probability(sequence, current_motif, j), sequence->background_site_probability[i][j]);

                    print_motif_model(motif_models[i]);
                }
            } else {
                sequence->site_probability_ratio[i][j] = 0;
            }
//            printf("%.20Lf ", sequence->site_probability_ratio[i][j]);
        }
//        printf("\n");
    }
//    printf("calculated site probability ratios\n");

    /**
     *  Calculate the site probabilities used for sampling.
     */
    for (j = 0; j < max_sites; j++) {
        sequence->site_probability[j][0] = 0;

        for (k = 1; k < sequence->length; k++) {
            sequence->site_probability[j][k] = sequence->site_probability[j][k - 1]; //this is the probability for the background of site j,k

            for (l = 0; l < number_motifs; l++) {
                current_motif = motif_models[l];

                motif_contribution = 0;
                if (possible_end_position(sequence, current_motif, k)) {
                    if (j == 0) { 
                        //previous probability is always 1 when j == 0
                        motif_contribution = sequence->site_probability_ratio[l][k];

                        sequence->motif_probability_contribution[j][k][l] = motif_contribution;
                        sequence->site_probability[j][k] += motif_contribution;
                    } else if (k - current_motif->motif_width < 0) {
                        //previous probability is 0
                        continue; 
                    } else {
                        motif_contribution = sequence->site_probability[j - 1][k - current_motif->motif_width] * sequence->site_probability_ratio[l][k]; 
                        sequence->motif_probability_contribution[j][k][l] = motif_contribution;
                        sequence->site_probability[j][k] += motif_contribution;
                    }
                }

                if (sequence->site_probability[j][k] < 0) {
                    fprintf(stderr, "ERROR: calculating site probabilities, probability became negative. file [%s] line [%d]\n", __FILE__, __LINE__);
                    fprintf(stderr, "sequence->site_probability[%d][%d]: %Lf, sequence->motif_probability_contribution[%d][%d][%d]: %Lf\n", j, k, sequence->site_probability[j][k], j, k, l, sequence->motif_probability_contribution[j][k][l]);
                    fprintf(stderr, "motif contribution: %Lf\n", motif_contribution);
                    if (j == 0) {
                        fprintf(stderr, "j == 0\n");
                    } else {
                        fprintf(stderr, "j == %d, sequence->site_probability[%d][%d]: %Lf, sequence->site_probability_ratio[%d][%d]: %Lf\n", j, j-1, k - current_motif->motif_width, sequence->site_probability[j-1][k-current_motif->motif_width], l, k, sequence->site_probability_ratio[l][k]);
                        fprintf(stderr, "site_probability_ratio[%d][%d]: %.30Lf\n", l, k-1, sequence->site_probability_ratio[l][k-1]);
                    }
                    exit(0);
                }
            }
        }
    }
}



void read_sequence_file(char *filename) {
    FILE *sequence_file;
    char current_line[200];
    int current_length, line_length, name_length;
    int i;
    Sequence *current_sequence;
    char c;

#ifdef _BOINC_    
    int retval;
    char input_path[512];

    retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
    if (retval) {

        fprintf(stderr, "ERROR, Could not resolve filename for reading sequences: [%s]\n", filename);
        fprintf(stderr, "ERROR on line [%d], file [%s]\n", __LINE__, __FILE__);
        exit(0);
    }

    sequence_file = boinc_fopen(input_path, "r");
#else
    sequence_file = fopen(filename, "r");
#endif
    if (sequence_file == NULL) {
        fprintf(stderr, "ERROR, could not open filename to read sequence file: [%s]\n", filename);
        fprintf(stderr, "ERROR on line [%d], file [%s]\n", __LINE__, __FILE__);
        exit(0);
    }

    number_sequences = 0;

//    fprintf(stderr, "mallocing sequences\n");
    sequences = (Sequence**)malloc(sizeof(Sequence*));
//    fprintf(stderr, "malloced sequences\n");


    current_sequence = NULL;
    current_length = 0;

    while (NULL != fgets(current_line, 200, sequence_file)) {

        if (strcmp(current_line, "") == 0 || strcmp(current_line, "\n") == 0) {
//            fprintf(stdout, "skipping line: [%s]\n", current_line);
            continue;

        } else if (current_line[0] == '>') {
            //start a new sequence

//            fprintf(stdout, "starting new sequence with line: [%s]\n", current_line);

            number_sequences++;
//            fprintf(stderr, "reallocing sequences\n");
            sequences = (Sequence**)realloc(sequences, number_sequences * sizeof(Sequence*));
//            fprintf(stderr, "realloced sequences\n");

            sequences[number_sequences - 1] = (Sequence*)malloc(sizeof(Sequence));
            current_sequence = sequences[number_sequences - 1];

            name_length = strlen(current_line) - 1;
            current_sequence->name = (char*)malloc(sizeof(char) * (name_length));

            for (i = 1; i < name_length; i++) current_sequence->name[i-1] = current_line[i];
            current_sequence->name[name_length-1] = '\0';

            current_sequence->nucleotides = NULL;
            current_length = 0;

//            fprintf(stdout, "reading new sequence, named: [%s]\n", current_sequence->name);

        } else {
            //append the current line to the sequence

//            fprintf(stdout, "appending line to sequence: [%s]\n", current_line);

            line_length = strlen(current_line);
//            printf("current line is     [%s], length [%d]\n", current_line, line_length);

            for (i = 0; i < line_length; i++) {
                c = toupper(current_line[i]);
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'X') {
                    current_line[i] = '\0';
                    line_length = i;
                    break;
                }
            }
//            printf("current line is     [%s], length [%d]\n", current_line, line_length);

//            current_line[line_length] = '\0';
            if (current_length == 0) {
//                fprintf(stderr, "mallocing nucleotides\n");
                current_sequence->nucleotides = (char*)malloc(line_length + 1);
//                fprintf(stderr, "malloced nucleotides\n");
            } else {
//                fprintf(stderr, "reallocing nucleotides\n");
                current_sequence->nucleotides = (char*)realloc(current_sequence->nucleotides, current_length + line_length + 1);
//                fprintf(stderr, "realloced nucleotides\n");
            }

            for (i = 0; i < line_length; i++) current_sequence->nucleotides[current_length + i] = current_line[i];
            current_sequence->nucleotides[current_length + line_length] = '\0';

            current_length += line_length;
            current_sequence->length = current_length;

            for (i = 0; i < current_sequence->length; i++) {
                c = toupper(current_sequence->nucleotides[i]);
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'X') {
                    printf("error in sequence, character at [%d] is not ACGTX: [%c][%d]\n", i, c, (int)c);
                    exit(0);
                }
            }
//            printf("nucleotides  is now [%s], length [%d]\n\n", current_sequence->nucleotides, current_sequence->length);
        }

    }

    fclose(sequence_file);

//    printf("read %d sequences\n", number_sequences);
//    for (i = 0; i < number_sequences; i++) {
//        printf("[%s]\n", sequences[i]->name);
//        printf("[%d] %s\n", sequences[i]->length, sequences[i]->nucleotides);
//    }

    /**
     *  Convert the sequences to uppercase to make things easier (so we don't have to check for lowercase as well).
     */
    sequences_to_uppercase();

    /**
     *  Calculate the background probabilities since we only need to do this once.
     */
    calculate_background_nucleotide_probabilities();

    /**
     *  Testing: print out the sequences
     */
    /*
    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];
        printf("name: [%s]\n", current_sequence->name);
        printf("nucleotides: [%s]\n", current_sequence->nucleotides);
        printf("length: [%d]\n", current_sequence->length);
        printf("strlen(nucleotides): [%d]\n", (int)strlen(current_sequence->nucleotides));

        printf("background_probability[A]: %.15lf\n", current_sequence->background_probability[A_POSITION]);
        printf("background_probability[C]: %.15lf\n", current_sequence->background_probability[C_POSITION]);
        printf("background_probability[G]: %.15lf\n", current_sequence->background_probability[G_POSITION]);
        printf("background_probability[T]: %.15lf\n", current_sequence->background_probability[T_POSITION]);

        if (strlen(current_sequence->nucleotides) != current_sequence->length) {
            //sanity check
            fprintf(stderr, "ERROR! strlen(nucleotides) != length!\n");
            exit(0);
        }

        printf("\n");
    }
    */
}

void write_sites_to_file(FILE* file, const char *delimiter) {
    int i, k;
    int already_printed;

    already_printed = 0;
    for (i = 0; i < number_sequences; i++) {
        for (k = 0; k < max_sites; k++) {
            if (sequences[i]->sampled_sites[k]->end_position >= 0) {
                if (already_printed) fprintf(file, ":");
                fprintf(file, "%d,%d", sequences[i]->sampled_sites[k]->motif_model, sequences[i]->sampled_sites[k]->end_position);
                already_printed = 1;
            }   
        }   
        fprintf(file, "%s", delimiter);
        already_printed = 0;
    }   
}

