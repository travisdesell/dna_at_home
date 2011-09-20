#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "structs.h"
#include "sequences.h"
#include "motif_models.h"
#include "sampling.h"

#include "../../mersenne_twister/dSFMT.h"


extern int max_sites;
extern int max_shift_distance;
extern int shift_possibilities;

extern int number_sequences;
extern Sequence **sequences;

extern int number_motifs;
extern MotifModel **motif_models;

MotifModel **shifted_motif_models;

void initialize_shifting() {
    int i;

    /**
     *  Note that the number of shifted motif models is not the number of motifs -- we're having a different model for each shift possibility
     */
    shifted_motif_models = (MotifModel**)malloc(sizeof(MotifModel*) * number_motifs);
    for (i = 0; i < number_motifs; i++) {
        shifted_motif_models[i] = (MotifModel*)malloc(sizeof(MotifModel));
        initialize_motif_model(shifted_motif_models[i], motif_models[i]->motif_width, motif_models[i]->type);
    }
}

void free_shifting() {
    int i;

    for (i = 0; i < number_motifs; i++) {
        free_motif_model(shifted_motif_models[i]);
        free(shifted_motif_models[i]);
    }
    free(shifted_motif_models);
}


void attempt_model_shift(int motif_model, int shift) {
    int i, j;
    Sequence *current_sequence;
    Sample *sample, *shifted_sample;
    long double ratio;
    long double rand;
    long double shifted_probability, unshifted_probability;
    int print_ratio;
    int temp_model_type;

    print_ratio = 0;

    zero_counts(shifted_motif_models[motif_model]);

//    printf("zeroed counts\n");

//    printf("\n\nshift: %d\n", shift);
    ratio = 0.0;

    for (i = 0; i < number_sequences; i++) {
//        printf("shifting sequence: %d\n", i);

        current_sequence = sequences[i];

        for (j = 0; j < max_sites; j++) {
            sample = current_sequence->sampled_sites[j];
            shifted_sample = current_sequence->shifted_sites[j];

            shifted_sample->motif_model = sample->motif_model;
            shifted_sample->end_position = sample->end_position;

//            printf("checking to see if shifted_sample->motif_model < 0: [%d]\n", shifted_sample->motif_model);

            if (shifted_sample->motif_model < 0) continue;
            if (shifted_sample->motif_model != motif_model) continue;

            shifted_sample->end_position += shift;

            if (!possible_end_position(current_sequence, motif_models[motif_model], shifted_sample->end_position)) {
//                printf("\nshift: %d, motif: %d\n", shift, motif_model);
//                printf("not shifting because sequence [%d] shifting [%d] from position [%d] to [%d] is not a possible end position\n", i, shift, sample->end_position, shifted_sample->end_position);
                shifted_sample->end_position -= shift;
                continue; //We can't shift this model by this much
            }

            /**
             *  Check to see if the shifted sample won't overlap another one.
             *  Since samples are taken in order, we know that the end positions are ascending and only need to check the adjacent samples
             */
            if (shift < 0 && j > 0 && current_sequence->sampled_sites[j-1]->motif_model > 0 && motif_model != current_sequence->sampled_sites[j-1]->motif_model && (shifted_sample->end_position - motif_models[motif_model]->motif_width) <= current_sequence->sampled_sites[j-1]->end_position) {
//                printf("\nshift: %d, motif: %d\n", shift, motif_model);
//                printf("not shifting because motif model [%d] != previous samples motif model [%d], and current samples end position - width [%d] <= previous samples end positon [%d]\n", motif_model, current_sequence->sampled_sites[j-1]->motif_model, shifted_sample->end_position - motif_width, current_sequence->sampled_sites[j-1]->end_position);
                shifted_sample->end_position -= shift;
                continue; // we cant shift this site to the left because a site would overlap another site (if it was the same model, that would be shfited by the same amount)

            } else if (shift > 0 && j < (max_sites - 1) && current_sequence->sampled_sites[j+1]->motif_model > 0 && motif_model != current_sequence->sampled_sites[j+1]->motif_model && shifted_sample->end_position > (current_sequence->sampled_sites[j+1]->end_position - motif_models[motif_model]->motif_width)) {
//                printf("\nshift: %d, motif: %d\n", shift, motif_model);
//                printf("not shifting because motif model [%d] != next samples motif model [%d], and current samples end position [%d] > next samples end positon - width [%d]\n", motif_model, current_sequence->sampled_sites[j+1]->motif_model, shifted_sample->end_position, current_sequence->sampled_sites[j+1]->end_position - motif_width);
                shifted_sample->end_position -= shift;
                continue; // we cant shift this site to the right because a site would overlap another site (if that was the same model, that would be shifted by the same amount)

            }
        }
    }

    //NEED TO UPDATE THE MOTIF MODEL
    for (i = 0; i < number_sequences; i++) {
        increment_model_counts(shifted_motif_models[motif_model], sequences[i], sequences[i]->shifted_sites, motif_model);
    }
    
    temp_model_type = shifted_motif_models[motif_model]->type;
    if (temp_model_type == MODEL_TYPE_FORWARD || temp_model_type == MODEL_TYPE_REVERSE) shifted_motif_models[motif_model]->type = MODEL_TYPE_NORMAL;

    update_motif_model(shifted_motif_models[motif_model]);

    shifted_motif_models[motif_model]->type = temp_model_type;


//    printf("\nSHIFT: %d\n", shift);
//    printf("SHIFTED:\n");
//    print_motif_model(shifted_motif_models[motif_model]);
//    printf("\nUNSHIFTED:\n");
//    print_motif_model(motif_models[motif_model]);

    shifted_probability = 0.0;
    unshifted_probability = 0.0;

    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];

        for (j = 0; j < max_sites; j++) {
            sample = current_sequence->sampled_sites[j];
            shifted_sample = current_sequence->shifted_sites[j];

            if (sample->motif_model == motif_model) {
//                printf("\n\nmotif_model: %d, sample->motif_model: %d, shifted_sample->motif_model: %d\n", motif_model, sample->motif_model, shifted_sample->motif_model);
//                printf("sample->end_position: %d, shifted_sample->end_position: %d\n", sample->end_position, shifted_sample->end_position);
//                printf("\nSHIFTED foreground_probability: %.20Lf, background_probability: %.20Lf\n", foreground_probability(current_sequence, shifted_motif_models[motif_model], shifted_sample->end_position), current_sequence->background_site_probability[motif_model][shifted_sample->end_position]);

                shifted_probability += log( foreground_probability(current_sequence, shifted_motif_models[motif_model], shifted_sample->end_position) / current_sequence->background_site_probability[motif_model][shifted_sample->end_position]);
                unshifted_probability += log( foreground_probability(current_sequence, motif_models[motif_model], sample->end_position) / current_sequence->background_site_probability[motif_model][sample->end_position]);
//                unshifted_probability += log ( current_sequence->site_probability_ratio[motif_model][sample->end_position] );
//                printf("incremented ratio to: %.20Lf\n", ratio);

//                printf("SHIFTED ratio: %.20Lf\n", foreground_probability(current_sequence, shifted_motif_models[motif_model], shifted_sample->end_position) / current_sequence->background_site_probability[motif_model][shifted_sample->end_position]);
//                printf("UNSHIFTED ratio: %.20Lf\n", current_sequence->site_probability_ratio[motif_model][sample->end_position]);

//                ratio -= log( current_sequence->site_probability_ratio[motif_model][sample->end_position] );
//                printf("ratio: %.20Lf\n", ratio);
            }
        }
    }

    ratio = shifted_probability - unshifted_probability;

    if (shifted_probability < unshifted_probability ) {
        /**
         *  In this case, the probability of taking the shift is < the probability of the current site in this case we 
         *  will move to the shift with probability = p(shift)/p(noshift), which is exp(ratio)
         */
        ratio = exp(ratio);
        rand = (long double)dsfmt_gv_genrand_open_close();
//        printf("moving to shift with probability: %.20Lf, rand: %.20Lf\n", ratio, rand);

        if (rand >= ratio) return; //don't make the shift
    }

    for (i = 0; i < number_sequences; i++) {
        current_sequence = sequences[i];

        for (j = 0; j < max_sites; j++) {
            current_sequence->sampled_sites[j]->end_position = current_sequence->shifted_sites[j]->end_position;
            current_sequence->sampled_sites[j]->motif_model = current_sequence->shifted_sites[j]->motif_model;
        }
    }

    copy_motif_model(shifted_motif_models[motif_model], motif_models[motif_model]);
}

void attempt_shifting() {
    int i, shift;

    shift = (int)(dsfmt_gv_genrand_open_close() * (double)max_shift_distance * 2.0);

    /**
     *  Select the shift randomly, but != 0
     */
//    printf("shift: %d\n", shift);
    shift -= max_shift_distance;
//    printf("shift - max_shift_distance: %d\n", shift);

    if (shift >= 0) shift++; 
//    printf("shift adjusted: %d\n", shift);

//    printf("shifting: %d\n", shift);

    for (i = 0; i < number_motifs; i++) {
//        printf("shifting model: %d\n", i);
        attempt_model_shift(i, shift);
    }
}
