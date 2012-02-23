#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "structs.hpp"
#include "sequences.hpp"
#include "motif_models.hpp"
#include "sampling.hpp"

#include "../../mersenne_twister/dSFMT.h"

using namespace std;

vector<MotifModel> shifted_motif_models;

void initialize_shifting(vector<MotifModel> &motif_models) {
    /**
     *  Note that the number of shifted motif models is not the number of motifs -- we're having a different model for each shift possibility
     */
    for (unsigned int i = 0; i < motif_models.size(); i++) {
        shifted_motif_models.push_back( MotifModel(motif_models.at(i).type, motif_models.at(i).motif_width) );
    }
}

void attempt_model_shift(int motif_model, int shift, vector<Sequence*> *sequences, vector<MotifModel> &motif_models, PhylogenyTree *phylogeny_tree) {
    Sample *sample, *shifted_sample, *adjacent_sample;
    long double shifted_probability, unshifted_probability;
    int temp_model_type;

    MotifModel *shifted_motif_model = &(shifted_motif_models.at(motif_model));

    shifted_motif_model->zero_counts();

    /**
     *  Need to shift in the opposite direction if the model is REVERSE
     */
    if (shifted_motif_model->type == MODEL_TYPE_REVERSE) shift = -shift;

//    printf("\n\nshift: %d\n", shift);

    bool made_shift = false;

    for (unsigned int i = 0; i < sequences->size(); i++) {
        Sequence *sequence = sequences->at(i);

        for (unsigned int j = 0; j < sequence->sampled_sites.size(); j++) {
            sample = &(sequence->sampled_sites.at(j));
            if (sample->motif_model < 0) continue;
            if (sample->motif_model != motif_model) continue;

            shifted_sample = &(sequence->shifted_sites.at(j));
            shifted_sample->motif_model = sample->motif_model;
            shifted_sample->end_position = sample->end_position;

            shifted_sample->end_position += shift;

            if (!sequence->possible_end_position(motif_models.at(motif_model), shifted_sample->end_position)) {
//                printf("\nshift: %d, motif: %d\n", shift, motif_model);
//                printf("not shifting because sequence [%d] shifting [%d] from position [%d] to [%d] is not a possible end position\n", i, shift, sample->end_position, shifted_sample->end_position);
                shifted_sample->end_position -= shift;
                continue; //We can't shift this model by this much
            }

            /**
             *  Check to see if the shifted sample won't overlap another one.
             *  Since samples are taken in order, we know that the end positions are ascending and only need to check the adjacent samples
             */
            if (shift < 0 && j > 0) {
                adjacent_sample = &(sequence->sampled_sites.at(j-1));

                if (adjacent_sample->motif_model > 0 && (shifted_sample->end_position - motif_models.at(motif_model).motif_width) <= adjacent_sample->end_position) {
    //                printf("\nshift: %d, motif: %d\n", shift, motif_model);
    //                printf("not shifting because motif model [%d] != previous samples motif model [%d], and current samples end position - width [%d] <= previous samples end positon [%d]\n", motif_model, sequence->sampled_sites[j-1]->motif_model, shifted_sample->end_position - motif_width, sequence->sampled_sites[j-1]->end_position);
                    shifted_sample->end_position -= shift;
                    continue; // we cant shift this site to the left because a site would overlap another site (if it was the same model, that would be shfited by the same amount)
                }

            } else if (shift > 0 && j < (sequence->sampled_sites.size() - 1)) {
                adjacent_sample = &(sequence->sampled_sites.at(j+1));

                if (adjacent_sample->motif_model > 0 && shifted_sample->end_position > (adjacent_sample->end_position - motif_models.at(motif_model).motif_width)) {
    //                printf("\nshift: %d, motif: %d\n", shift, motif_model);
    //                printf("not shifting because motif model [%d] != next samples motif model [%d], and current samples end position [%d] > next samples end positon - width [%d]\n", motif_model, sequence->sampled_sites[j+1]->motif_model, shifted_sample->end_position, sequence->sampled_sites[j+1]->end_position - motif_width);
                    shifted_sample->end_position -= shift;
                    continue; // we cant shift this site to the right because a site would overlap another site (if that was the same model, that would be shifted by the same amount)
                }
            }
            made_shift = true;
        }
    }

    if (!made_shift) {
//        cerr << "no changes made to sequence, was not possible to shift any samples." << endl;
        return;
    }
    //NEED TO UPDATE THE MOTIF MODEL
    for (unsigned int i = 0; i < sequences->size(); i++) {
        shifted_motif_model->increment_model_counts(sequences->at(i), sequences->at(i)->shifted_sites, motif_model);
    }
    
    temp_model_type = shifted_motif_model->type;
    if (temp_model_type == MODEL_TYPE_FORWARD || temp_model_type == MODEL_TYPE_REVERSE) shifted_motif_model->type = MODEL_TYPE_NORMAL;

    shifted_motif_model->update_motif_model();
    shifted_motif_model->type = temp_model_type;

//    printf("\nSHIFT: %d\n", shift);
//    printf("SHIFTED:\n");
//    shifted_motif_models[motif_model].print(cerr);
//    printf("\nUNSHIFTED:\n");
//    motif_models[motif_model].print(cerr);

    shifted_probability = 0.0;
    unshifted_probability = 0.0;

    for (unsigned int i = 0; i < sequences->size(); i++) {
        Sequence *sequence = sequences->at(i);

        for (unsigned int j = 0; j < sequence->sampled_sites.size(); j++) {
            sample = &(sequence->sampled_sites[j]);
            shifted_sample = &(sequence->shifted_sites[j]);

            if (sample->motif_model == motif_model) {
//                printf("\n\nmotif_model: %d, sample->motif_model: %d, shifted_sample->motif_model: %d\n", motif_model, sample->motif_model, shifted_sample->motif_model);
//                printf("sample->end_position: %d, shifted_sample->end_position: %d\n", sample->end_position, shifted_sample->end_position);
//                printf("\nSHIFTED foreground_probability: %.20Lf, background_probability: %.20Lf\n", foreground_probability(sequence, shifted_motif_models[motif_model], shifted_sample->end_position), sequence->background_site_probability[motif_model][shifted_sample->end_position]);

                if (phylogeny_tree == NULL) {
                    shifted_probability += log( sequence->foreground_probability(*shifted_motif_model, shifted_sample->end_position) / sequence->background_site_probability.at(motif_model).at(shifted_sample->end_position));
                    unshifted_probability += log( sequence->foreground_probability(motif_models.at(motif_model), sample->end_position) / sequence->background_site_probability.at(motif_model).at(sample->end_position));
                } else {
                    shifted_probability += log( sequence->foreground_probability_phylogeny(*shifted_motif_model, shifted_sample->end_position, phylogeny_tree) );
                    unshifted_probability += log( sequence->foreground_probability_phylogeny(motif_models.at(motif_model), sample->end_position, phylogeny_tree) );
                }
//                unshifted_probability += log ( sequence->site_probability_ratio[motif_model][sample->end_position] );
//                printf("incremented ratio to: %.20Lf\n", ratio);

//                printf("SHIFTED ratio: %.20Lf\n", foreground_probability(sequence, shifted_motif_models[motif_model], shifted_sample->end_position) / sequence->background_site_probability[motif_model][shifted_sample->end_position]);
//                printf("UNSHIFTED ratio: %.20Lf\n", sequence->site_probability_ratio[motif_model][sample->end_position]);

//                ratio -= log( sequence->site_probability_ratio[motif_model][sample->end_position] );
//                printf("ratio: %.20Lf\n", ratio);
            }
        }
    }

    long double ratio = shifted_probability - unshifted_probability;
    long double rand = 0;

    if (shifted_probability < unshifted_probability ) {
        /**
         *  In this case, the probability of taking the shift is < the probability of the current site in this case we 
         *  will move to the shift with probability = p(shift)/p(noshift), which is exp(ratio)
         */
        ratio = exp(ratio);
        rand = (long double)dsfmt_gv_genrand_open_close();

        if (rand >= ratio) return; //don't make the shift
    }
//  printf("moved to shift with probability: %.20Lf, rand: %.20Lf\n", ratio, rand);

//    cout << "\nPERFORMING SHIFT of " << shift << "!\n" << endl;

    for (unsigned int i = 0; i < sequences->size(); i++) {
        Sequence *sequence = sequences->at(i);

        for (unsigned int j = 0; j < sequence->sampled_sites.size(); j++) {
            if (sequence->shifted_sites.at(j).motif_model == motif_model) {
  //              cout << endl << "shifted sequence->sampled_sites[j].end_position to: " << sequence->sampled_sites.at(j).end_position << " from " << sequence->shifted_sites.at(j).end_position << endl;
 //               cout << "shifted sequence->sampled_sites[j].motif_model to: " << sequence->sampled_sites.at(j).motif_model << " from " << sequence->shifted_sites.at(j).motif_model << endl;
 //               cout << "sequence->sampled_sites.size(): " << sequence->sampled_sites.size() << endl;
 //               cout << "sequence->shifted_sites.size(): " << sequence->shifted_sites.size() << endl;
                sequence->sampled_sites.at(j).end_position = sequence->shifted_sites.at(j).end_position;
                sequence->sampled_sites.at(j).motif_model = sequence->shifted_sites.at(j).motif_model;
            }
        }
    }

//    cout << "copying motif models" << endl;

//    cout << "shifted motif model:" << endl;
//    shifted_motif_model->print(cout);
//    cout << "motif model [" << motif_model << "]:" << endl;
//    motif_models.at(motif_model).print(cout);

    copy_motif_model(*shifted_motif_model, motif_models.at(motif_model));
    /**
     *  If a shift happens, we need to zero out the model counts and recalculate them.
     */
}

void attempt_shifting(int max_shift_distance, vector<Sequence*> *sequences, vector<MotifModel> &motif_models, PhylogenyTree *phylogeny_tree) {
    int shift = (int)(dsfmt_gv_genrand_open_close() * (double)max_shift_distance * 2.0);

    /**
     *  Select the shift randomly, but != 0
     */
//    printf("shift: %d\n", shift);
    shift -= max_shift_distance;
//    printf("shift - max_shift_distance: %d\n", shift);

    if (shift >= 0) shift++; 
//    printf("shift adjusted: %d\n", shift);

//    printf("shifting: %d\n", shift);

    for (unsigned int i = 0; i < motif_models.size(); i++) {
        attempt_model_shift(i, shift, sequences, motif_models, phylogeny_tree);
    }

    /**
     *  Should call attempt model shift on the model itself
     */

}
