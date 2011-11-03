#ifndef SEQUENCES_H
#define SEQUENCES_H

#include <iostream>
#include <vector>

typedef class Sample;
typedef class Sequence;

#include "structs.hpp"
#include "motif_models.hpp"

class Sample {
    public:
        int end_position;
        int motif_model;

        Sample() {
            end_position = -1;
            motif_model = -1;
        }

        Sample(int _motif_model, int _end_position) {
            motif_model = _motif_model;
            end_position = _end_position;
        }
};

class Sequence {
    public:
        string name;

        string nucleotides;
        vector<int> distance_to_previous_invalid;      // used for calculating if a position is a possible end position for a motif model

        vector< vector<long double> > possible_site_counts;  // counts the number of possible ways to place [1 .. max_sites] motif model sites within the sequence

        vector<long double> possible_sites;                //the values of the last column of possible site counts
        long double total_possible_sites;           //the sum of possible_sites

        vector< vector<long double> > background_site_probability;
        vector< vector<long double> > site_probability_ratio;

        vector< vector<long double> > site_probability;
        vector< vector< vector<long double> > > motif_probability_contribution;

        vector<Sample> sampled_sites;
        vector<Sample> shifted_sites;

        vector< vector<long> > accumulated_samples;

        Sequence(string sequence_information, string nucleotides, int max_sites, int number_motifs);

        int possible_end_position(MotifModel &motif_model, unsigned int position);

        long double background_probability(double *background_nucleotide_probability, MotifModel &motif_model, int end_position);
        long double foreground_probability(MotifModel &motif_model, int end_position);

        void calculate_distance_to_invalid();
        void calculate_background_site_probabilities(vector<MotifModel> &motif_models);
        void calculate_site_probabilities(vector<MotifModel> &motif_models, int max_sites);

        void zero_accumulated_samples();

        void sample_uniform_random(vector<MotifModel> &motif_models);
        void sample_uniform_random_helper(int number_sites, int position, vector<MotifModel> &motif_models);
        void resample_from_models_helper(int numer_sites, int position, vector<MotifModel> &motif_models);
        void resample_from_models(vector<MotifModel> &motif_models);
};

void read_sequences(vector<Sequence> &sequences, string sequence_filename, int max_sites, int number_motifs);
void calculate_background_nucleotide_probabilities(vector<Sequence> &sequences);
void calculate_background_site_probabilities(vector<Sequence> &sequence);

void calculate_site_counts(vector<Sequence> &sequences, vector<MotifModel> &motif_models, vector<double> &blocks, int max_sites);

void read_sequence_file(char *filename);

#endif
