#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>

#include "phylogeny.hpp"
#include "structs.hpp"
#include "util.hpp"

#ifdef _USE_GSL_
    #include "./lee_phylogeny/SubstitutionMatrix.h"
#endif

using namespace std;

#ifdef PHYLOGENY_TEST
int main(int argc, char**argv) {
    vector<string> arguments(argv, argv + argc);

    /**
     *  Lee's code uses A T C G
     *  BG: (A: 0.20, T: 0.20, C: 0.30, G: 0.30)
     *  FG: (A: .05, T: 0.10, C: 0.80, G: 0.05)
     */
    vector<double> *bg_equilibrium = new vector<double>(4, 0.0);
    bg_equilibrium->at(0) = 0.2;    
    bg_equilibrium->at(1) = 0.3;
    bg_equilibrium->at(2) = 0.3;
    bg_equilibrium->at(3) = 0.2;

    vector<double> *fg_equilibrium = new vector<double>(4, 0.0);
    fg_equilibrium->at(0) = 0.05;
    fg_equilibrium->at(1) = 0.80;
    fg_equilibrium->at(2) = 0.05;
    fg_equilibrium->at(3) = 0.10;

    double tt_factor = 3.0;

    int max_sites = 1;
    int max_shift_distance = 0;
    int motif_info_size = 1;

    /**
     * read sequences from file: ../tests/phylogeny_test/simple.seq
     * read phylogeny tree from file: ../tests/phylogeny_test/simple.phylip
     */
    vector<Sequence*> *sequences = new vector<Sequence*>(0, NULL);
    read_sequences(sequences, "../tests/phylogeny_test/simple.seq", max_sites, motif_info_size /*motif_info.size() is the number of motifs*/, max_shift_distance);

    cout << "****************************" << endl;
    cout << "TEST WITH 2 SEQUENCES" << endl;
    cout << "****************************" << endl;
    cout << endl;

    ifstream in1("../tests/phylogeny_test/simple.phylip");
    PhylogenyTree *phylogeny_tree = new PhylogenyTree(in1, sequences);

    long double bg_prob = phylogeny_tree->get_background_probability(bg_equilibrium, tt_factor, 0);
    cout << endl;
    cout << "bg_prob: " << bg_prob << endl;

    long double fg_prob = phylogeny_tree->get_foreground_probability(fg_equilibrium, 0);

    cout << endl;
    cout << "simple:" << phylogeny_tree->get_phylip_string() << endl;
    cout << "bg_prob: " << bg_prob << ", fg_prob: " << fg_prob << endl;

//    bg_prob = phylogeny_tree->get_background_probability(bg_equilibrium, tt_factor, 0);
//    fg_prob = phylogeny_tree->get_foreground_probability(fg_equilibrium, 0);
//    cout << "RECALCULATION TEST: bg_prob: " << bg_prob << ", fg_prob: " << fg_prob << endl;

    /**
     * read sequences from file: ../tests/phylogeny_test/simple2.seq
     * read phylogeny tree from file: ../tests/phylogeny_test/simple2.phylip
     */
    for (unsigned int i = 0; i < sequences->size(); i++) delete sequences->at(i);
    delete sequences;
    sequences = new vector<Sequence*>(0, NULL);
    read_sequences(sequences, "../tests/phylogeny_test/simple2.seq", max_sites, motif_info_size /*motif_info.size() is the number of motifs*/, max_shift_distance);

    delete phylogeny_tree;

    cout << endl;
    cout << endl;
    cout << "****************************" << endl;
    cout << "TEST WITH 3 SEQUENCES" << endl;
    cout << "****************************" << endl;
    cout << endl;

    ifstream in2("../tests/phylogeny_test/simple2.phylip");
    phylogeny_tree = new PhylogenyTree(in2, sequences);

    bg_prob = phylogeny_tree->get_background_probability(bg_equilibrium, tt_factor, 0);
    cout << endl;
    cout << "   bg_prob: " << bg_prob << endl;

    fg_prob = phylogeny_tree->get_foreground_probability(fg_equilibrium, 0);

    cout << endl;
    cout << "simple2:" << phylogeny_tree->get_phylip_string() << endl;
    cout << "bg_prob: " << bg_prob << ", fg_prob: " << fg_prob << endl;

//    bg_prob = phylogeny_tree->get_background_probability(bg_equilibrium, tt_factor, 0);
//    fg_prob = phylogeny_tree->get_foreground_probability(fg_equilibrium, 0);
//    cout << "RECALCULATION TEST: bg_prob: " << bg_prob << ", fg_prob: " << fg_prob << endl;

    for (unsigned int i = 0; i < sequences->size(); i++) delete sequences->at(i);
    delete sequences;
    delete phylogeny_tree;
}

#endif

/**
 *  Helper functions (for parsing the phylip file)
 */

bool is_whitespace(char c) {
    return c == ' ' || c == '\n' || c == '\r' || c == '\t';
}

void skip_whitespace(string str, int &pos) {
    while (pos < (int)str.size() && is_whitespace( str.at(pos) )) {
        pos++;
    }
}

string make_whitespace(int length) {
    string s(length, ' ');
    return s;
}

int get_next_colon(string str, int pos) {
    while (pos < (int)str.size() && str.at(pos) != ':') {
        pos++;
    }
    return pos;
}

/**
 *  PhylogenyTreeNode methods
 */

PhylogenyTreeNode::~PhylogenyTreeNode() {
    for (unsigned int i = 0; i < children->size(); i++) {
        delete children->at(i);
    }
    delete nucleotide_probabilities;
    delete children;
    leaf_sequence = NULL;
}

PhylogenyTreeNode::PhylogenyTreeNode(string phylip_string, int &start, vector<Sequence*> *sequences, vector<bool> &sequence_match) {
    nucleotide_probabilities = new vector<long double>(4, 0.0);
    children = new vector<PhylogenyTreeNode*>();
    leaf_sequence = NULL;

//    cout << "parsing node at:" << endl;
//    cout << "   " << phylip_string << endl;
//    cout << "   " << make_whitespace(start) << "^" << endl;

    if (phylip_string.at(start) == '(') {   //this node has children
        leaf = false;

        while (phylip_string.at(start) != ')') {
            start++;
            skip_whitespace(phylip_string, start);

            children->push_back( new PhylogenyTreeNode(phylip_string, start, sequences, sequence_match) );

            if (phylip_string.at(start) != ',' && phylip_string.at(start) != ')') {
                cerr << "ERROR: incorrectly formatted phylip string: " << endl;
                cerr << "   expected ',' or ')' after node: " << start << endl;
                cerr << "   " << phylip_string << endl;
                cerr << "   " << make_whitespace(start) << "^" << endl;
                cerr << "   occurred in [" << __FILE__ << "] on line: " << __LINE__ << endl;
                exit(1);
            }
        }
        start++;
        skip_whitespace(phylip_string, start);

        if (phylip_string.at(start) == ';') {
            evolution_time = 0;
            sequence_name = "ROOT";
            return;
        } else if (phylip_string.at(start) != ':') {
            cerr << "ERROR: incorrectly formatted phylip string: " << endl;
            cerr << "   no colon found for child starting at index: " << start << endl;
            cerr << "   " << phylip_string << endl;
            cerr << "   " << make_whitespace(start) << "^" << endl;
            cerr << "   occurred in [" << __FILE__ << "] on line: " << __LINE__ << endl;
            exit(1);
        } 
        sequence_name = "INTERNAL";

    } else {    //this is a leaf node
        leaf = true;
        int colon_pos = phylip_string.find(':', start);

        if (colon_pos > (int)phylip_string.size()) {
                cerr << "ERROR: incorrectly formatted phylip string: " << endl;
                cerr << "   expected colon after: " << start << endl;
                cerr << "   " << phylip_string << endl;
                cerr << "   " << make_whitespace(start) << "^" << endl;
                cerr << "   occurred in [" << __FILE__ << "] on line: " << __LINE__ << endl;
                exit(1);
        }
        sequence_name = phylip_string.substr(start, colon_pos - start);

        for (unsigned int i = 0; i < sequences->size(); i++) {
            int space_pos = sequences->at(i)->name.find(' ');
            string tmp_name = sequences->at(i)->name.substr(1, space_pos - 1);

//            cout << " comparing sequence names: '" << sequence_name << "' to '" << tmp_name << "'" << endl;
            if (tmp_name.compare(sequence_name) == 0) {
                leaf_sequence = sequences->at(i);
                sequence_match.at(i) = true;
                break;
            }
        }

        if (leaf_sequence == NULL) {
            cerr << "ERROR: could not find matching sequence for phylip tree leaf: " << sequence_name << endl;
            cerr << "   occurred in [" << __FILE__ << "] on line: " << __LINE__ << endl;
            exit(1);
        }
        start = colon_pos;
    }
    start++;
    skip_whitespace(phylip_string, start);

    int comma_pos = phylip_string.find(',', start);
    int close_paren_pos = phylip_string.find(')', start);

    int end = (comma_pos > start  && comma_pos < close_paren_pos) ? comma_pos : close_paren_pos;

    if (end > (int)phylip_string.size()) {
            cerr << "ERROR: incorrectly formatted phylip string: " << endl;
            cerr << "   no closing ',' or ')' found after index: " << start << endl;
            cerr << "   " << phylip_string << endl;
            cerr << "   " << make_whitespace(start) << "^" << endl;
            cerr << "   occurred in [" << __FILE__ << "] on line: " << __LINE__ << endl;
            exit(1);
    }

    evolution_time = atof( phylip_string.substr(start + 1, (end - start) - 1).c_str() );

//    cout << "sequence_name: " << sequence_name << endl;
//    cout << "evolution_time: " << evolution_time << endl;
    start = end;
}

void PhylogenyTreeNode::print(int indent) {
    cout << make_whitespace(indent) << "sequence_name: " << sequence_name << endl;
    cout << make_whitespace(indent) << "evolution_time: " << evolution_time << endl;
    cout << make_whitespace(indent) << "nucleotide_probabilities: " << vector_to_string<long double>(nucleotide_probabilities) << endl;
    cout << endl;

    for (unsigned int i = 0; i < children->size(); i++) {
        children->at(i)->print(indent+1);
    }
}

void PhylogenyTreeNode::set_leaf_nucleotide_probabilities(int sequence_position) {
    if (is_leaf()) {    //Lee's code uses A T C G
        nucleotide_probabilities->at(0) = 0.0;
        nucleotide_probabilities->at(1) = 0.0;
        nucleotide_probabilities->at(2) = 0.0;
        nucleotide_probabilities->at(3) = 0.0;

        if (leaf_sequence->nucleotides.at(sequence_position) == 'A') nucleotide_probabilities->at(0) = 1;
        if (leaf_sequence->nucleotides.at(sequence_position) == 'T') nucleotide_probabilities->at(1) = 1;
        if (leaf_sequence->nucleotides.at(sequence_position) == 'C') nucleotide_probabilities->at(2) = 1;
        if (leaf_sequence->nucleotides.at(sequence_position) == 'G') nucleotide_probabilities->at(3) = 1;
    } else {
        for (unsigned int i = 0; i < children->size(); i++) {
            children->at(i)->set_leaf_nucleotide_probabilities(sequence_position);
        }
    }
}

#ifdef _USE_GSL_
void PhylogenyTreeNode::calculate_nucleotide_probabilities(const Wadsworth::SubstitutionMatrix &substitution_matrix) {
    if (!is_leaf()) {
        vector< Wadsworth::SubstitutionMatrix::SMMatrix > *responses = new vector< Wadsworth::SubstitutionMatrix::SMMatrix>(children->size(), Wadsworth::SubstitutionMatrix::SMMatrix(ALPHABET_LENGTH, ALPHABET_LENGTH, Wadsworth::MATRIX_UNINITIALIZED));

        for (unsigned int i = 0; i < children->size(); i++) {
            children->at(i)->calculate_nucleotide_probabilities(substitution_matrix);
            substitution_matrix.getSubstitutionMatrix(children->at(i)->get_evolution_time(), responses->at(i));
        }

        long double current_probability;
        for (unsigned int i = 0; i < ALPHABET_LENGTH; i++) {
            nucleotide_probabilities->at(i) = 1.0;
            for (unsigned int j = 0; j < children->size(); j++) {
                current_probability = 0.0;
                for (unsigned int k = 0; k < ALPHABET_LENGTH; k++) {
                    current_probability += responses->at(j)(i,k) * children->at(j)->nucleotide_probabilities->at(k);
                }
                //current_probability should be (c[j]->prob dot response[j][i])
                nucleotide_probabilities->at(i) *= current_probability;
            }
        }

        delete responses;
    }
}
#endif

/**
 *  PhylogenyTree methods
 */

PhylogenyTree::~PhylogenyTree() {
    delete root;
}


PhylogenyTree::PhylogenyTree(string _phylip_string, vector<Sequence*> *sequences) {
    init(_phylip_string, sequences);
}

PhylogenyTree::PhylogenyTree(ifstream &in, vector<Sequence*> *sequences) {
    string _phylip_string;
    string line;

    while (in.good()) {
        getline(in, line);
        _phylip_string.append(line);
    }
    init(_phylip_string, sequences);
}

void PhylogenyTree::init(string _phylip_string, vector<Sequence*> *sequences) {
    phylip_string = _phylip_string;

    int position = 0;
    skip_whitespace(phylip_string, position);

    if (phylip_string.at(position) != '(') {
        cerr << "ERROR: incorrectly formatted phylip string: " << endl;
        cerr << "   does not begin with open parenthesis on index " << position << endl;
        cerr << "   " << phylip_string << endl;
        cerr << "   " << make_whitespace(position) << "^" << endl;
        cerr << "   occurred in [" << __FILE__ << "] on line: " << __LINE__ << endl;
        exit(1);
    }

    vector<bool> sequence_match(sequences->size(), false);
//    cout << "\nparsed phylogeny tree:" << endl << endl;

    root = new PhylogenyTreeNode(phylip_string, position, sequences, sequence_match);

    bool error = false;
    for (unsigned int i = 0; i < sequences->size(); i++) {
        if (!sequence_match.at(i)) {
            cerr << "ERROR: unmatched sequence in phylip tree: " << sequences->at(i)->name << endl;
            error = true;
        }
    }
    if (error) {
        cerr << "   occurred in [" << __FILE__ << "] on line: " << __LINE__ << endl;
        exit(1);
    }

//    print();
//    cout << "." << endl << endl;
}

void PhylogenyTree::set_background_equilibrium(vector<double> *background_equilibrium) {
#ifdef _USE_GSL_
    Wadsworth::SubstitutionMatrix::SMColVector bg_equilibrium(ALPHABET_LENGTH, 1, 0.0); //Lee's code uses A T C G
    bg_equilibrium(0,0) = background_equilibrium->at(A_POSITION);
    bg_equilibrium(1,0) = background_equilibrium->at(T_POSITION);
    bg_equilibrium(2,0) = background_equilibrium->at(C_POSITION);
    bg_equilibrium(3,0) = background_equilibrium->at(G_POSITION);

//    cout << "set background equilibrium: [" << bg_equilibrium(0,0) << ", " << bg_equilibrium(1,0) << ", " << bg_equilibrium(2,0) << ", " << bg_equilibrium(3,0) << "]" << endl;

    substitution_matrix.setEquilibriumAndTTFactor( bg_equilibrium, tt_factor );
#endif
}

void PhylogenyTree::set_foreground_equilibrium(vector<double> *foreground_equilibrium) {
#ifdef _USE_GSL_
    Wadsworth::SubstitutionMatrix::SMColVector fg_equilibrium(ALPHABET_LENGTH, 1, 0.0); //Lee's code uses A T C G
    fg_equilibrium(0,0) = foreground_equilibrium->at(A_POSITION);
    fg_equilibrium(1,0) = foreground_equilibrium->at(T_POSITION);
    fg_equilibrium(2,0) = foreground_equilibrium->at(C_POSITION);
    fg_equilibrium(3,0) = foreground_equilibrium->at(G_POSITION);

//    cout << "foreground equilibrium: [" << fg_equilibrium(0,0) << ", " << fg_equilibrium(1,0) << ", " << fg_equilibrium(2,0) << ", " << fg_equilibrium(3,0) << "]" << endl;

    substitution_matrix.setHB98Foreground(fg_equilibrium);
#endif
}

void PhylogenyTree::set_leaf_nucleotide_probabilties(int sequence_position) {
    root->set_leaf_nucleotide_probabilities(sequence_position);
}

long double PhylogenyTree::get_background_probability(vector<double> *background_equilibrium, int sequence_position) {
    set_background_equilibrium(background_equilibrium);

    set_leaf_nucleotide_probabilties(sequence_position);
//    cout << "after set nucleotide probabilities:" << endl;
//    print();
//    cout << "." << endl;

#ifdef _USE_GSL_
    root->calculate_nucleotide_probabilities(substitution_matrix);
#endif
//    cout << "after calculate nucleotide probabilities:" << endl;
//    print();
//    cout << "." << endl;

    long double bg_probability = 0;
    bg_probability += root->get_nucleotide_probability(0) * background_equilibrium->at(A_POSITION);
    bg_probability += root->get_nucleotide_probability(1) * background_equilibrium->at(T_POSITION);
    bg_probability += root->get_nucleotide_probability(2) * background_equilibrium->at(C_POSITION);
    bg_probability += root->get_nucleotide_probability(3) * background_equilibrium->at(G_POSITION);

    return bg_probability;
}

long double PhylogenyTree::get_foreground_probability(vector<double> *foreground_equilibrium) {
    set_foreground_equilibrium(foreground_equilibrium);

#ifdef _USE_GSL_
    root->calculate_nucleotide_probabilities(substitution_matrix);
#endif
//    cout << "after calculate nucleotide probabilities:" << endl;
//    print();

    double fg_probability = 0;
    fg_probability += root->get_nucleotide_probability(0) * foreground_equilibrium->at(A_POSITION);
    fg_probability += root->get_nucleotide_probability(1) * foreground_equilibrium->at(T_POSITION);
    fg_probability += root->get_nucleotide_probability(2) * foreground_equilibrium->at(C_POSITION);
    fg_probability += root->get_nucleotide_probability(3) * foreground_equilibrium->at(G_POSITION);

//    cout << "fg_prob: " << fg_probability << endl;

    return fg_probability;
}

void PhylogenyTree::set_tt_factor(double _tt_factor) {
    tt_factor = _tt_factor;
}

void PhylogenyTree::print() {
    root->print(0);
}
