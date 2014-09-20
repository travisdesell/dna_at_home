#ifndef PHYLOGENY_H
#define PHYLOGENY_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

typedef class PhylogenyTree;

#include "sequences.hpp"

#ifdef _USE_GSL_
    #include "./lee_phylogeny/SubstitutionMatrix.h"
#endif

using namespace std;

class PhylogenyTreeNode {
    protected:
        /**
         *  When calculating the probability a new nucleotide, these probabilities are set (at the leaf nodes) to, e.g. 0 0 1 0 -- if the nucleotide is C (given an A T C G order)
         */
        vector<long double> *nucleotide_probabilities;

        double evolution_time;

        vector<PhylogenyTreeNode*> *children;

        Sequence* leaf_sequence;

        int sequence_id;
        string sequence_name;
        bool leaf;

    public:
        inline double get_evolution_time()      { return evolution_time; }
        inline int get_sequence_id()            { return sequence_id; }
        inline string get_sequence_name()          { return sequence_name; }
        inline bool is_leaf()                   { return leaf; }
        inline bool is_internal()               { return !leaf; }

        ~PhylogenyTreeNode();

        PhylogenyTreeNode(string phylip_string, int &start, vector<Sequence*> *sequences, vector<bool> &sequence_match);

        long double get_nucleotide_probability(int i) { return nucleotide_probabilities->at(i); }

        void set_leaf_nucleotide_probabilities(int sequence_position);
#ifdef _USE_GSL_
        void calculate_nucleotide_probabilities(const Wadsworth::SubstitutionMatrix &substitution_matrix);
#endif

        void print(int indent);
};

class PhylogenyTree {
    protected:
        string phylip_string;
        PhylogenyTreeNode *root;
        double tt_factor;

#ifdef _USE_GSL_
        Wadsworth::SubstitutionMatrix substitution_matrix;
#endif

    public:
        string get_phylip_string() { return phylip_string; }

        PhylogenyTree(string _phylip_string, vector<Sequence*> *sequences);
        PhylogenyTree(ifstream &in, vector<Sequence*> *sequences);
        ~PhylogenyTree();

        void set_tt_factor(double _tt_factor);

        void init (string _phylip_string, vector<Sequence*> *sequences);

        void set_background_equilibrium(vector<double> *background_equilibrium);
        void set_foreground_equilibrium(vector<double> *foreground_equilibrium);

        void set_leaf_nucleotide_probabilties(int sequence_position);

        long double get_background_probability(vector<double> *background_equilibrium, int sequence_position);
        long double get_foreground_probability(vector<double> *foreground_equilibrium);

        void print();
};

#endif
