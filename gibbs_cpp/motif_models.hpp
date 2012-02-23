#ifndef MOTIF_MODELS_H
#define MOTIF_MODELS_H

typedef class MotifModel;

#include "structs.hpp"
#include "sequences.hpp"

#define MODEL_TYPE_NORMAL 0
#define MODEL_TYPE_REVERSE 1
#define MODEL_TYPE_FORWARD 2
#define MODEL_TYPE_PALINDROMIC 3

const string MODEL_TYPES[4] = { "NORMAL", "REVERSE", "FORWARD", "PALINDROMIC" };

using namespace std;

class MotifModel {
    public: 
        int type;
        int motif_width;

        /**
         * counts is sized [motif_width][ALPHABET_LENGTH], it contains the sum of each letter found by the motif model at its sampled positions
         */
        vector< vector<long> > counts;

        /**
         * nucleotide_probabilities is sized [motif_width][ALPHABET_LENGTH], it contains the probability of finding a letter according to counts
         */
        vector< vector<double> > nucleotide_probabilities;
        vector< vector<double> > reverse_nucleotide_probabilities;

        MotifModel(int _type, int _motif_width);
        MotifModel(string motif_info);
        MotifModel(int _type, vector< vector<double> > &_nucleotide_probabilities);

        ~MotifModel() {}

        void increment_counts_for_sample(Sequence *sequence, const Sample &sample);
        void increment_model_counts(Sequence *sequence, const vector<Sample> &samples, int motif_model_number);

        void zero_counts();

        void print(ostream &out_stream);
        void print_short(ostream &out_stream);

        void update_motif_model();

        string generate_possible_nucleotides();
};                                                                                            

void update_motif_models(vector<MotifModel> &motif_models);
void update_motif_model_reverse_complement(MotifModel &forward_motif_model, MotifModel &reverse_motif_model);

void initialize_motif_models(vector<MotifModel> &motif_models, vector<string> &motif_info);

void copy_motif_model(const MotifModel &source, MotifModel &destination);

void increment_counts(vector<MotifModel> &motif_models, Sequence *sequence);
void decrement_counts(vector<MotifModel> &motif_models, Sequence *sequence);

void read_motifs_from_file(string filename, vector<MotifModel> &motif_models);

#endif
