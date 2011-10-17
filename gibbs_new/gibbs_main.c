#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"
#include "sequences.h"
#include "motif_models.h"
#include "sampling.h"
#include "time.h"
#include "shifting.h"
#include "checkpoint.h"

#ifdef _WIN32
    #include "boinc_win.h"
    #include "str_util.h"
#endif

#ifdef _BOINC_
    #include "diagnostics.h"
    #include "util.h"
    #include "filesys.h"
    #include "boinc_api.h"
    #include "mfile.h"
#endif

#include "../../mersenne_twister/dSFMT.h"


#define SITES_CHECKPOINT_FILE "gibbs_sites_checkpoint.txt"
#define SAMPLES_CHECKPOINT_FILE "gibbs_samples_checkpoint.txt"
#define SAMPLES_OUTPUT_FILE "accumulated_samples.txt"

int max_sites;

int number_sequences;
int max_shift_distance;
Sequence **sequences;

double background_nucleotide_probability[ALPHABET_LENGTH];

double foreground_pseudocounts = 0.28;
double background_pseudocounts = 5.0;


int number_motifs;
MotifModel **motif_models;

/**
 *  accumulated_samples is sized [number_sequences][number_motifs][sequence_length[i]]
 */
int ***accumulated_samples;


int main(int number_arguments, char** arguments) {
    int i, j, k, retval, iteration;
    int burn_in_period, sample_period;
    int shift_period;
    int seed;
    int sites_from_arguments;
    Sequence *current_sequence;
    Sample **current_sample;
    double *blocks;
    double count_percentage;
    int print_current_sites, print_accumulated_samples, print_current_sites_frequency;
    double print_best_sites;
    double progress;
    int total_number_nucleotides;

    retval = 0;
    #ifdef _BOINC_
        #ifdef BOINC_APP_GRAPHICS
            #if defined(_WIN32) || defined(__APPLE)
                retval = boinc_init_graphics(worker);
            #else
                retval = boinc_init_graphics(worker, argv[0]);
            #endif
        #else
            retval = boinc_init();
        #endif

        if (retval) exit(retval);
    #endif

    seed = -1;
    max_shift_distance = -1;
    shift_period = -1;
    print_best_sites = 0.0;
    print_current_sites = 0;
    print_current_sites_frequency = 0;
    print_accumulated_samples = 0;
    max_sites = -1;
    iteration = 0;
    sample_period = 0;
    burn_in_period = 0;
    sites_from_arguments = 0;

    number_motifs = 0;

    for (i = 0; i < number_arguments; i++) {
        fprintf(stderr, "argument [%d]: [%s]\n", i, arguments[i]);

        if (!strcmp(arguments[i], "--motifs")) {
            number_motifs = atoi(arguments[++i]);
            i++;
            i += initialize_motif_models(i, arguments);

        } else if (!strcmp(arguments[i], "--sequence_file")) {
            read_sequence_file(arguments[++i]);

        } else if (!strcmp(arguments[i], "--burn_in_period")) {
            burn_in_period = atoi(arguments[++i]);

        } else if (!strcmp(arguments[i], "--sample_period")) {
            sample_period = atoi(arguments[++i]);

        } else if (!strcmp(arguments[i], "--max_sites")) {
            max_sites = atoi(arguments[++i]);
            blocks = (double*)malloc(sizeof(double) * (max_sites + 1));
            for (j = 0; j < max_sites + 1; j++) blocks[j] = 1.0 / ((double)max_sites + 1.0);

        } else if (!strcmp(arguments[i], "--seed")) {
            seed = atoi(arguments[++i]);

        } else if (!strcmp(arguments[i], "--enable_shifting")) {
            max_shift_distance = atoi(arguments[++i]);
            shift_period = atoi(arguments[++i]);

        } else if (!strcmp(arguments[i], "--print_current_sites")) {
            print_current_sites = 1;

        } else if (!strcmp(arguments[i], "--print_current_sites_frequency")) {
            print_current_sites_frequency = atoi(arguments[++i]);

        } else if (!strcmp(arguments[i], "--print_accumulated_samples")) {
            print_accumulated_samples = 1;

        } else if (!strcmp(arguments[i], "--print_best_sites")) {
            print_best_sites = atof(arguments[++i]);

        }

    }

    initialize_sequences_data();

    blocks = NULL;
    for (i = 0; i < number_arguments; i++) {
//        fprintf(stderr, "argument [%d]: [%s]\n", i, arguments[i]);

        if (!strcmp(arguments[i], "--blocks")) {
            if (max_sites < 0) {
                printf("ERROR, max sites not specified.\n\n");
//                print_usage();
                exit(0);
            }

            blocks = (double*)malloc(sizeof(double) * (max_sites + 1));
            for (j = 0; j < max_sites + 1; j++) blocks[j] = atof(arguments[++i]);

        } else if (!strcmp(arguments[i], "--current_sites")) {
            sites_from_arguments = 1;
            read_sites(arguments[++i]);

        }
    }

    if (blocks == NULL) {
        fprintf(stderr, "ERROR, blocks == NULL!\n");
        fflush(stderr);
        exit(0);
    }

    fprintf(stderr, "blocks:");
    for (i = 0; i < max_sites + 1; i++) fprintf(stderr, " %lf", blocks[i]);
    fprintf(stderr, "\n");

    if (seed < 0) seed = time(NULL);


    if (shift_period > 0) initialize_shifting();

    calculate_distance_to_invalid();
    calculate_site_counts(blocks);

    /**
     * select initial random samples
     */

    accumulated_samples = (int***)malloc(sizeof(int**) * number_sequences);
    for (i = 0; i < number_sequences; i++) {
        accumulated_samples[i] = (int**)malloc(sizeof(int*) * number_motifs);

        for (j = 0; j < number_motifs; j++) {
            accumulated_samples[i][j] = (int*)malloc(sizeof(int) * sequences[i]->length);

            for (k = 0; k < sequences[i]->length; k++) {
                accumulated_samples[i][j][k] = 0;
            }
        }
    }

    if (!read_sites_from_checkpoint(SITES_CHECKPOINT_FILE, &seed, &iteration)) {
        fprintf(stderr, "seeding: %d\n", seed);
        dsfmt_gv_init_gen_rand(seed);

        //Sample from the sequences uniformly at random if there is no checkpoint file.
        if (!sites_from_arguments) {
            for (i = 0; i < number_sequences; i++) {
                sample_uniform_random(sequences[i]);
            }
        } else {
            fprintf(stderr, "sites were from arguments\n");
        }

    } else {
        if (iteration >= burn_in_period) {
            read_accumulated_samples(SAMPLES_CHECKPOINT_FILE);
        }

        fprintf(stderr, "seeding: %d\n", seed + iteration);
        dsfmt_gv_init_gen_rand(seed + iteration);
    }

    for (i = 0; i < number_sequences; i++) {
        calculate_background_site_probabilities(sequences[i]);
    }

    total_number_nucleotides = 0;
    for (i = 0; i < number_sequences; i++) {
        total_number_nucleotides += sequences[i]->length;
    }
    fprintf(stderr, "%d sequences with %d total base pairs.\n", number_sequences, total_number_nucleotides);

    /**
     * calculate the counts for all the samples -- we don't need to do this at 
     * every step, we can subtract the left out sequence from the counts then 
     * add the newly sampled counts
     */
//    fprintf(stderr, "incrementing counts from checkpoint for [%d] sequences:\n", number_sequences);
    for (i = 0; i < number_sequences; i++) {
//        fprintf(stderr, "incrementing counts for sequence [%d][%s][%s]\n", i, sequences[i]->name, sequences[i]->nucleotides);

//        for (j = 0; j < max_sites; j++) {
//            fprintf(stderr, "\tsample[%d]->motif_model: %d, sample->[%d]->end_position: %d\n", j, sequences[i]->sampled_sites[j]->motif_model, j, sequences[i]->sampled_sites[j]->end_position);
//        }

        increment_all_counts(sequences[i]);

//        fprintf(stderr, "\n");
    }
    fprintf(stderr, "incremented counts from checkpoint for [%d] sequences.\n", number_sequences);

    printf("burn in period: %d, sample period: %d\n", burn_in_period, sample_period);
    /**
     * do the burn in walk
     */

    if (print_current_sites_frequency > 0) {
        fprintf(stderr, "<current_sites>\n");
        fprintf(stderr,"<iteration>0</iteration>\n");
        write_sites_to_file(stderr, ".\n");
        fprintf(stderr, "</current_sites>\n");
    }

    for (i = iteration; i < burn_in_period + sample_period; i++) {

        if (i > 0 && shift_period > 0 && (i % shift_period) == 0) {
            attempt_shifting();
        }

        for (j = 0; j < number_sequences; j++) {
            current_sequence = sequences[j];
            current_sample = current_sequence->sampled_sites;

            /**
             * leave one sequence out and re-calculate the motif models
             * we can do this by decrementing the counts of the left-out sequence, and recalculating the models
             */
            decrement_counts(current_sequence);
           
            /**
             * update the models
             */
            for (k = 0; k < number_motifs; k++) {
                switch (motif_models[k]->type) {
                    case MODEL_TYPE_FORWARD:
                        update_motif_model_reverse_complement(motif_models[k], motif_models[k + 1]);
                        k++;
                        break;
                    default:
                        update_motif_model(motif_models[k]);
                        break;
                }
            }

            /**
             * resample within the sequence left out
             */
            calculate_site_probabilities(current_sequence);

            resample_from_models(current_sequence);
            
            /**
             * update the counts with the new sample
             */
            increment_all_counts(current_sequence);

            for (k = 0; k < number_motifs; k++) {
                switch (motif_models[k]->type) {
                    case MODEL_TYPE_FORWARD:
                        update_motif_model_reverse_complement(motif_models[k], motif_models[k + 1]);
                        k++;
                        break;
                    default:
                        update_motif_model(motif_models[k]);
                        break;
                }
            }

            if (i >= burn_in_period) {
                //add the samples taken this iteration to the saved samples

                for (k = 0; k < max_sites; k++) {
                    if (current_sample[k]->end_position < 0 || current_sample[k]->motif_model < 0) continue;
                    accumulated_samples[j][current_sample[k]->motif_model][current_sample[k]->end_position]++;
                }
            }
        
            progress = ((double)i)/((double)(burn_in_period + sample_period));
            progress += (((double)j)/((double)number_sequences) * (1.0/((double)(burn_in_period + sample_period))));

#ifdef _BOINC_
            boinc_fraction_done(progress);
#else
            fprintf(stdout, "\r%lf", progress);
#endif
        }

//#ifdef _BOINC_
//        if (boinc_time_to_checkpoint()) {
//            write_sites(SITES_CHECKPOINT_FILE, seed, i + 1);
//            if (i >= burn_in_period) write_accumulated_samples(SAMPLES_CHECKPOINT_FILE);

//            dsfmt_gv_init_gen_rand(seed + i + 1);

//            boinc_checkpoint_completed();
//        }
//#else
        if (i % 5000 == 0 && i != 0) {
            write_sites(SITES_CHECKPOINT_FILE, seed, i + 1);
            
            if (i >= burn_in_period) write_accumulated_samples(SAMPLES_CHECKPOINT_FILE);

            dsfmt_gv_init_gen_rand(seed + i + 1);
#ifdef _BOINC_
            boinc_checkpoint_completed();
#endif
        }

        if (print_current_sites_frequency > 0 && (((i+1) % print_current_sites_frequency) == 0)) {
            fprintf(stderr, "<current_sites>\n");
            fprintf(stderr,"<iteration>%d</iteration>\n", i + 1);
            write_sites_to_file(stderr, ".\n");
            fprintf(stderr, "</current_sites>\n");
        }

//#endif
    }

    printf("\n");
//        printf("samples for sequence [%s]\n", sequences[i]->nucleotides);

    if (print_best_sites > 0) {
        for (j = 0; j < number_motifs; j++) {
            for (i = 0; i < number_sequences; i++) {

                for (k = 0; k < sequences[i]->length; k++) {
                    count_percentage = ((double)accumulated_samples[i][j][k]) / ((double)sample_period);

                    if (count_percentage > print_best_sites) {
                        printf("%5d,%2d%8d ", i, j, (k + 1) - motif_models[j]->motif_width);
                        print_sample_and_nearest(sequences[i], motif_models[j], k, 5);
                        printf(" %7d %2.2lf %s\n", k, count_percentage, sequences[i]->name);
                    }
                }
            }
            printf("\n");
        }
    }

    if (print_current_sites) {
        fprintf(stderr, "<current_sites>\n");
        write_sites_to_file(stderr, ".\n");
        fprintf(stderr, "</current_sites>\n");
    }

    if (print_accumulated_samples) {
        fprintf(stderr, "<current_samples>\n");
        write_accumulated_samples_to_file(stderr);
        fprintf(stderr, "</current_samples>\n");
    }

    if (sample_period > 0) {
        write_accumulated_samples(SAMPLES_OUTPUT_FILE);
    }

    if (shift_period > 0) free_shifting();

    for (i = 0; i < number_motifs; i++) {
        print_motif_model(motif_models[i]);
    }

    free_motif_models();

    free_sequences();

    for (i = 0; i < number_sequences; i++) {
        for (j = 0; j < number_motifs; j++) {
            free(accumulated_samples[i][j]);
        }
        free(accumulated_samples[i]);
    }
    free(accumulated_samples);

    #ifdef _BOINC_
        boinc_finish(0);
    #endif

    return 0;
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode){
    LPSTR command_line;
    char* argv[100];
    int argc;

    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}
#endif
