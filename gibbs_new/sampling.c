#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "sampling.h"
#include "motif_models.h"
#include "sequences.h"

#include "../../mersenne_twister/dSFMT.h"

extern int max_sites;

extern int number_motifs;
extern MotifModel **motif_models;

//extern int number_sequences;
//extern Sequence **sequences;

void sample_uniform_random_helper(Sequence *sequence, int number_sites, int position) {
    int i;
    long double rand;
    long double total;
    long double no_site_probability;
    long double site_probability;
    MotifModel *motif_model;
    Sample *sample;

    if (number_sites < 0) return; //we've sampled all the site we're going to
    if (position - 1 <= 0) return;  //we've reached the beginning of the sequence

//    printf("\nSAMPLE UNIFORM RANDOM HELPER, position: [%d], number_sites: [%d]\n", position, number_sites);

    total = sequence->possible_site_counts[number_sites][position];
    no_site_probability = sequence->possible_site_counts[number_sites][position - 1];

    rand = dsfmt_gv_genrand_open_close() * total;

//    printf("rand: [%Lf], no_site_probability: [%Lf], total: [%Lf], position: [%d]\n", rand, no_site_probability, total, position);

    if (rand < no_site_probability) {
        sample_uniform_random_helper(sequence, number_sites, position - 1);
        return;

    } else {
        rand -= no_site_probability;


        for (i = 0; i < number_motifs; i++) {
            motif_model = motif_models[i];
//            printf("motif width: %d\n", motif_model->width);
//
//            printf("sampling for motif: [%d], number sites: [%d]\n", i, number_sites);

            if (number_sites == 0) {
                site_probability = 1;
            } else {
                site_probability = sequence->possible_site_counts[number_sites - 1][position - motif_model->motif_width];
            }

//            printf("rand: %Lf, total: %Lf, no_site_probability: %Lf, site_probability: %Lf, number sites: %d, position: %d\n", rand, total, no_site_probability, site_probability, number_sites, position);

//            printf("rand: [%Lf], site_probability[%d]: [%Lf], total: [%Lf]\n", rand, i, site_probability, total);

            if (rand < site_probability) {
//                printf("TAKING A SAMPLE\n");
                //take a sample
                sample = sequence->sampled_sites[number_sites];

                sample->end_position = position;
                sample->motif_model = i;

//                printf("calling helper with number sites: %d, position: %d\n", number_sites - 1, position - motif_width);
                sample_uniform_random_helper(sequence, number_sites - 1, position - motif_model->motif_width);
                return;
            } else {
                rand -= site_probability;
            }

        }

        fprintf(stderr, "\n\nERROR: sample uniform random helper reached a spot it should not have.  Tried all sites and no site probabilities and did not make a move");
        fprintf(stderr, "In file [%s], line [%d]\n", __FILE__, __LINE__);
        fprintf(stderr, "Sequence [%s]\n", sequence->nucleotides);
        fprintf(stderr, "position [%d], number sites [%d]\n", position, number_sites);
        exit(0);
    }
}

void sample_uniform_random(Sequence *sequence) {
    int i;
    long double rand;
    int number_sites;
    Sample *sample;

    for (i = 0; i < max_sites; i++) {
        sample = sequence->sampled_sites[i];

        sample->end_position = -1;
        sample->motif_model= -1;
    }

    //determine the number of sites to sample
    rand = (long double)dsfmt_gv_genrand_open_close() * (long double)sequence->total_possible_sites;

//    printf("total possible sites: %Lf\n", sequence->total_possible_sites);

    number_sites = 0;
    for (i = 0; i < max_sites + 1; i++) {
//        printf("sequence->possible_sites[%d]: %Lf\n", i, sequence->possible_sites[i]);

        if (rand < sequence->possible_sites[i]) {
            number_sites = i;
            break;
        } else {
            rand -= sequence->possible_sites[i];
        }
    }
//    printf("sampling %d sites\n", number_sites);

    if (number_sites == 0) return; // we aren't sampling any sites

    sample_uniform_random_helper(sequence, number_sites - 1, sequence->length - 1);
}




void sample_from_models_helper(Sequence *sequence, int number_sites, int position) {
    int i, j;
    long double rand;
    long double total;
    long double no_site_probability;
    long double site_probability;
    MotifModel *motif_model;
    Sample *sample;

    if (number_sites < 0) return; //we've sampled all the site we're going to
    if (position - 1 <= 0) return;  //we've reached the beginning of the sequence

//    printf("\nSAMPLE FROM MODEL HELPER, position: [%d], number_sites: [%d]\n", position, number_sites);

    total = sequence->site_probability[number_sites][position];
    no_site_probability = sequence->site_probability[number_sites][position - 1];

    rand = dsfmt_gv_genrand_open_close() * total;

//    printf("rand: [%.30Lf], number_sites: [%d], no_site_probability: [%.30Lf], total: [%.30Lf]\n", rand, number_sites, no_site_probability, total);

    if (rand < no_site_probability) {
        sample_from_models_helper(sequence, number_sites, position - 1);
        return;

    } else {
        rand -= no_site_probability;

        for (i = 0; i < number_motifs; i++) {
            motif_model = motif_models[i];
//            printf("motif width: %d\n", motif_model->width);

            site_probability = sequence->motif_probability_contribution[number_sites][position][i];

//            printf("rand: [%.15Lf], number_sites: [%d], site_probability[%d]: [%.15Lf], total: [%.15Lf]\n", rand, number_sites, i, site_probability, total);

            if (rand < site_probability) {
                //take a sample
                
//                printf("SAMPLED AT POSITION: %d\n", position);
                sample = sequence->sampled_sites[number_sites];

                sample->end_position = position;
                sample->motif_model = i;

                sample_from_models_helper(sequence, number_sites - 1, position - motif_model->motif_width);
                return;
            } else {
                rand -= site_probability;
            }
        }

        fprintf(stderr, "\n\nERROR: sample from models helper reached a spot it should not have.  Tried all sites and no site probabilities and did not make a move.\n");
        fprintf(stderr, "In file [%s], line [%d]\n", __FILE__, __LINE__);

        for (i = 0; i < number_motifs; i++) {
            print_motif_model(motif_models[i]);
            fprintf(stderr, "\n");
        }

        printf("\n\n");
        printf("Sequence:            [%s]\n", sequence->nucleotides);
    
        printf("Possible End Positions:\n");
        for (i = 0; i < number_motifs; i++) {
            printf("motif [%d], width [%d]: ", i, motif_models[i]->motif_width);
            for (j = 0; j < sequence->length; j++) {
                printf("%d", possible_end_position(sequence, motif_models[i], j));
            }
            printf("\n");
        }

        printf("distance to previous invalids:");
        for (i = 0; i < sequence->length; i++) {
            printf("%d ", sequence->distance_to_previous_invalid[i]);
        }
        printf("\n");

        printf("Counts:\n");
        for (i = 0; i < max_sites; i++) {
            for (j = 0; j < sequence->length; j++) {
                printf("%Lf ", sequence->possible_site_counts[i][j]);
            }
            printf("\n");
        }

        printf("Background Site Probabilities:\n");
        for (i = 0; i < number_motifs; i++) {
            for (j = 0; j < sequence->length; j++) {
                printf("%.20Lf ", sequence->background_site_probability[i][j]);
            }
            printf("\n");
        }

        printf("Site Probabilities:\n");
        for (i = 0; i < number_motifs; i++) {
            for (j = 0; j < sequence->length; j++) {
                printf("%Lf ", sequence->site_probability[i][j]);
            }
            printf("\n");
        }


        printf("position [%d], number sites [%d]\n", position, number_sites);
        exit(0);
    }
}



void resample_from_models(Sequence *sequence) {
    int i;
    long double rand;
    int number_sites;
    Sample *sample;

    /**
     *  Reset the samples
     */
    for (i = 0; i < max_sites; i++) {
        sample = sequence->sampled_sites[i];

        sample->end_position = -1;
        sample->motif_model= -1;
    }


    //determine the number of sites to sample
    rand = (long double)dsfmt_gv_genrand_open_close() * (long double)sequence->total_possible_sites;

    number_sites = 0;
    for (i = 0; i < max_sites + 1; i++) {
        if (rand < sequence->possible_sites[i]) {
            number_sites = i;
            break;
        } else {
            rand -= sequence->possible_sites[i];
        }
    }

    if (i == 0) return;    //if i == 0 we aren't sampling any sites

    sample_from_models_helper(sequence, number_sites - 1, sequence->length - 1);
}

void print_sample(Sequence *sequence, MotifModel *motif_model, int end_position) {
    int i;
    printf("[");
    for (i = 0; i < motif_model->motif_width; i++) {
        printf("%c", sequence->nucleotides[end_position - motif_model->motif_width + 1 + i]);
    }
    printf("]");
}

void print_sample_and_nearest(Sequence *sequence, MotifModel *motif_model, int end_position, int nearest) {
    int i, position;

    for (i = 0; i < nearest; i++) {
        position = end_position - motif_model->motif_width + 1 - nearest + i;
        if (position < 0) printf(" ");
        else printf("%c", tolower(sequence->nucleotides[position]));
    }
    printf(" ");

    for (i = 0; i < motif_model->motif_width; i++) {
        printf("%c", sequence->nucleotides[end_position - motif_model->motif_width + 1 + i]);
    }

    printf(" ");
    for (i = 0; i < nearest; i++) {
        position = end_position + 1 + i;
        if (position >= sequence->length) printf(" ");
        else printf("%c", tolower(sequence->nucleotides[position]));
    }
}


