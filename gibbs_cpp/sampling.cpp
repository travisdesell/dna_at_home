#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "structs.hpp"
#include "sampling.hpp"
#include "motif_models.hpp"
#include "sequences.hpp"

#include "../../mersenne_twister/dSFMT.h"

using namespace std;

void Sequence::sample_uniform_random_helper(int number_sites, int position, vector<MotifModel> &motif_models) {
    long double rand;
    long double total;
    long double no_site_prob;
    long double site_prob;

    if (number_sites < 0) return; //we've sampled all the site we're going to
    if (position - 1 <= 0) return;  //we've reached the beginning of the sequence

//    printf("\nSAMPLE UNIFORM RANDOM HELPER, position: [%d], number_sites: [%d]\n", position, number_sites);

    total = possible_site_counts.at(number_sites).at(position);
    no_site_prob = possible_site_counts.at(number_sites).at(position - 1);

    rand = dsfmt_gv_genrand_open_close() * total;

//    printf("rand: [%Lf], no_site_prob: [%Lf], total: [%Lf], position: [%d]\n", rand, no_site_prob, total, position);

    if (rand < no_site_prob) {
        sample_uniform_random_helper(number_sites, position - 1, motif_models);
        return;

    } else {
        rand -= no_site_prob;

        for (unsigned int i = 0; i < motif_models.size(); i++) {
            MotifModel *motif_model = &(motif_models.at(i));
//            printf("motif width: %d\n", motif_model->width);
//
//            printf("sampling for motif: [%d], number sites: [%d]\n", i, number_sites);

            if (number_sites == 0) {
                site_prob = 1;
            } else {
                site_prob = possible_site_counts.at(number_sites - 1).at(position - motif_model->motif_width);
            }

//            printf("rand: %Lf, total: %Lf, no_site_prob: %Lf, site_prob: %Lf, number sites: %d, position: %d\n", rand, total, no_site_prob, site_prob, number_sites, position);

//            printf("rand: [%Lf], site_prob[%d]: [%Lf], total: [%Lf]\n", rand, i, site_prob, total);

            if (rand < site_prob) {
//                printf("TAKING A SAMPLE\n");
                //take a sample
                Sample *sample = &(sampled_sites.at(number_sites));

                sample->end_position = position;
                sample->motif_model = i;

//                printf("calling helper with number sites: %d, position: %d\n", number_sites - 1, position - motif_width);
                sample_uniform_random_helper(number_sites - 1, position - motif_model->motif_width, motif_models);
                return;
            } else {
                rand -= site_prob;
            }

        }

        cerr << endl << endl << "ERROR: sample uniform random helper reached a spot it should not have.  Tried all sites and no site probabilities and did not make a move" << endl;
        cerr << "In file [" << __FILE__ << "], line [" << __LINE__ << "]" << endl;
        cerr << "Sequence [" << nucleotides << "]" << endl;
        cerr << "position [" << position << "], number sites [" << number_sites << "]" << endl;
        exit(0);
    }
}

void Sequence::sample_uniform_random(vector<MotifModel> &motif_models) {
    int number_sites;
    Sample *sample;

    for (unsigned int i = 0; i < sampled_sites.size(); i++) {
        sample = &(sampled_sites.at(i));

        sample->end_position = -1;
        sample->motif_model= -1;
    }

    //determine the number of sites to sample
    long double rand = (long double)dsfmt_gv_genrand_open_close() * (long double)total_possible_sites;

//    printf("total possible sites: %Lf\n", sequence->total_possible_sites);

    number_sites = 0;
    for (unsigned int i = 0; i < possible_sites.size(); i++) {
//        printf("sequence->possible_sites[%d]: %Lf\n", i, sequence->possible_sites[i]);

        if (rand < possible_sites.at(i)) {
            number_sites = i;
            break;
        } else {
            rand -= possible_sites.at(i);
        }
    }
//    printf("sampling %d sites\n", number_sites);

    if (number_sites == 0) return; // we aren't sampling any sites

    sample_uniform_random_helper(number_sites - 1, nucleotides.size() - 1, motif_models);
}

void Sequence::resample_from_models_helper(int number_sites, int position, vector<MotifModel> &motif_models) {
    unsigned int i, j;
    long double rand;
    long double current_site_probability;

    if (number_sites < 0) return; //we've sampled all the site we're going to
    if (position - 1 <= 0) return;  //we've reached the beginning of the sequence

//    printf("\nSAMPLE FROM MODEL HELPER, position: [%d], number_sites: [%d]\n", position, number_sites);

    long double total = site_probability.at(number_sites).at(position);
    long double no_site_probability = site_probability.at(number_sites).at(position - 1);

    rand = dsfmt_gv_genrand_open_close() * total;

//    printf("rand: [%.30Lf], number_sites: [%d], no_site_probability: [%.30Lf], total: [%.30Lf]\n", rand, number_sites, no_site_probability, total);

    if (rand < no_site_probability) {
        resample_from_models_helper(number_sites, position - 1, motif_models);
        return;

    } else {
        rand -= no_site_probability;

        for (i = 0; i < motif_models.size(); i++) {
            MotifModel *motif_model = &(motif_models.at(i));
//            printf("motif width: %d\n", motif_model->width);

            current_site_probability = motif_probability_contribution.at(number_sites).at(position).at(i);

//            printf("rand: [%.15Lf], number_sites: [%d], site_probability[%d]: [%.15Lf], total: [%.15Lf]\n", rand, number_sites, i, site_probability, total);

            if (rand < current_site_probability) {
                //take a sample
                
//                printf("SAMPLED AT POSITION: %d\n", position);
                Sample *sample = &(sampled_sites.at(number_sites));

                sample->end_position = position;
                sample->motif_model = i;

                resample_from_models_helper(number_sites - 1, position - motif_model->motif_width, motif_models);
                return;
            } else {
                rand -= current_site_probability;
            }
        }

        cerr << "ERROR: sample from models helper reached a spot it should not have.  Tried all sites and no site probabilities and did not make a move." << endl;
        cerr << "In file [" << __FILE__ << "], line [" << __LINE__ << "]" << endl;

        for (i = 0; i < motif_models.size(); i++) {
            motif_models.at(i).print(cerr);
            cerr << endl;
        }

        cerr << endl << endl;
        cerr << "Sequence:            [" << nucleotides << "]" << endl;
    
        cerr << "Possible End Positions:" << endl;
        for (i = 0; i < motif_models.size(); i++) {
            cerr << "motif [" << i << "], width [" << motif_models.at(i).motif_width << "]" ;
            for (j = 0; j < nucleotides.size(); j++) {
                cerr << possible_end_position(motif_models.at(i), j);
            }
            cerr << endl;
        }

        cerr << "distance to previous invalids:";
        for (i = 0; i < nucleotides.size(); i++) {
            cerr << distance_to_previous_invalid.at(i);
        }
        cerr << endl;

        cerr << "Counts:" << endl;
        for (i = 0; i < possible_site_counts.size(); i++) {
            for (j = 0; j < possible_site_counts.at(i).size(); j++) {
                cerr << possible_site_counts.at(i).at(j);
            }
            cerr << endl;
        }

        cerr << "Background Site Probabilities:" << endl;
        for (i = 0; i < background_site_probability.size(); i++) {
            for (j = 0; j < background_site_probability.at(i).size(); j++) {
                cerr << background_site_probability.at(i).at(j);
            }
            cerr << endl;
        }

        cerr << "Site Probabilities:";
        for (i = 0; i < site_probability.size(); i++) {
            for (j = 0; j < site_probability.at(i).size(); j++) {
                cerr << site_probability.at(i).at(j);
            }
            cerr << endl;
        }

        cerr << "position [" << position << "], number sites [" << number_sites << "]" << endl;
        exit(0);
    }
}

void Sequence::resample_from_models(vector<MotifModel> &motif_models) {
    unsigned int i;
    long double rand;
    int number_sites;
    Sample *sample;

    /**
     *  Reset the samples
     */
    for (i = 0; i < sampled_sites.size(); i++) {
        sample = &(sampled_sites.at(i));

        sample->end_position = -1;
        sample->motif_model= -1;
    }


    //determine the number of sites to sample
    rand = (long double)dsfmt_gv_genrand_open_close() * (long double)total_possible_sites;

    number_sites = 0;
    for (i = 0; i < possible_sites.size(); i++) {
        if (rand < possible_sites.at(i)) {
            number_sites = i;
            break;
        } else {
            rand -= possible_sites.at(i);
        }
    }

    if (i == 0) return;    //if i == 0 we aren't sampling any sites

    resample_from_models_helper(number_sites - 1, nucleotides.size() - 1, motif_models);
}

void print_sample(Sequence *sequence, MotifModel &motif_model, int end_position) {
    int i;
    printf("[");
    for (i = 0; i < motif_model.motif_width; i++) {
        printf("%c", sequence->nucleotides.at(end_position - motif_model.motif_width + 1 + i));
    }
    printf("]");
}

void print_sample_and_nearest(Sequence *sequence, MotifModel &motif_model, int end_position, int nearest) {
    int i, position;

    for (i = 0; i < nearest; i++) {
        position = end_position - motif_model.motif_width + 1 - nearest + i;
        if (position < 0) printf(" ");
        else printf("%c", tolower(sequence->nucleotides[position]));
    }
    printf(" ");

    for (i = 0; i < motif_model.motif_width; i++) {
        printf("%c", sequence->nucleotides[end_position - motif_model.motif_width + 1 + i]);
    }

    printf(" ");
    for (i = 0; i < nearest; i++) {
        position = end_position + 1 + i;
        if (position >= (int)sequence->nucleotides.size()) printf(" ");
        else printf("%c", tolower(sequence->nucleotides[position]));
    }
}


