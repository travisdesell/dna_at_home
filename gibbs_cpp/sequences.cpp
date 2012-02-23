#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#ifdef _BOINC_
    #include "diagnostics.h"
    #include "util.h"
    #include "filesys.h"
    #include "boinc_api.h"
    #include "mfile.h"
#endif

#include "motif_models.hpp"

#include "sequences.hpp"
#include "util.hpp"

using namespace std;

double background_nucleotide_probability[ALPHABET_LENGTH];

double foreground_pseudocounts = 0.28;
double background_pseudocounts = 5.0;

int get_nucleotide_position(char nucleotide) {
    switch (nucleotide) {
        case 'A': return A_POSITION;
        case 'C': return C_POSITION;
        case 'G': return G_POSITION;
        case 'T': return T_POSITION;
        default: return -1;
    }
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

Sequence::Sequence(string sequence_information, string _nucleotides, int max_sites, int number_motifs, int max_shift_distance) {
    name = sequence_information;
    nucleotides = _nucleotides;
    distance_to_previous_invalid = vector<int>((int)nucleotides.size(), 0);
    calculate_distance_to_invalid();

    possible_sites = vector<long double>(max_sites + 1, (long double)0.0);
    possible_site_counts = vector< vector<long double> >(max_sites, vector<long double>(nucleotides.size(), (long double)0.0));

    background_site_probability = vector< vector<long double> >(number_motifs, vector<long double>(nucleotides.size(), (long double)0.0));
    site_probability_ratio = vector< vector<long double> >(number_motifs, vector<long double>(nucleotides.size(), (long double)0.0));

    site_probability = vector< vector<long double> >(max_sites, vector<long double>(nucleotides.size(), (long double)0.0));
    motif_probability_contribution = vector< vector< vector<long double> > >(max_sites, vector< vector<long double> >(nucleotides.size(), vector<long double>(number_motifs, (long double)0.0)));

    sampled_sites = vector<Sample>(max_sites, Sample());
    if (max_shift_distance > 0) shifted_sites = vector<Sample>(max_sites, Sample());

    accumulated_samples = vector< vector<long> >(number_motifs, vector<long>(nucleotides.size(), 0));
}

int Sequence::possible_end_position(MotifModel &motif_model, unsigned int position) {
    if (position < 1 || position >= nucleotides.size() || distance_to_previous_invalid.at(position) < motif_model.motif_width) return 0;
    return 1;
}

void Sequence::calculate_distance_to_invalid() {
    if (is_valid_nucleotide(nucleotides.at(0))) {
        distance_to_previous_invalid.at(0) = 1;
    } else {
        distance_to_previous_invalid.at(0) = 0;
    }

    for (unsigned int i = 1; i < distance_to_previous_invalid.size(); i++) {
        if (is_valid_nucleotide(nucleotides.at(i))) {
            distance_to_previous_invalid.at(i) = distance_to_previous_invalid.at(i-1) + 1;
        } else {
            distance_to_previous_invalid.at(i) = 0;
        }
    }
}

void remove_carriage_return(string &s) {
    for (string::iterator it = s.begin(); it < s.end(); it++) {
        if (*it == '\r') s.erase(it);
    }
}

void read_sequences(vector<Sequence*> *sequences, string sequence_filename, int max_sites, int number_motifs, int max_shift_distance) {
    ifstream sequence_file(sequence_filename.c_str());
    if (sequence_file.is_open()) {
        string current_line;

        getline(sequence_file, current_line);
        remove_carriage_return(current_line);
        while (sequence_file.good() && !sequence_file.eof()) {
            /**
             *  File format is (for each sequence):
             *  >sequence information
             *  NUCLEOTIDES
             *  NUCLEOTIDES
             *  NUCLEOTIDES
             *  <blank line>
             */
            string sequence_information = current_line;

            string nucleotides;
            getline(sequence_file, current_line);
            remove_carriage_return(current_line);
//            cout << "next sequence line is: " << current_line << endl;
            do {
                nucleotides += current_line;
                getline(sequence_file, current_line);
                remove_carriage_return(current_line);

//                cout << "next sequence line is: " << current_line << endl;
            } while (sequence_file.good() && !sequence_file.eof() && current_line.compare("") != 0);
//            cout << "finished reading nucleotides." << endl;

            for (unsigned int i = 0; i < nucleotides.size(); i++) {
                /**
                 *  Convert the sequence to all uppercase
                 */
                nucleotides.at(i) = toupper( (unsigned char)nucleotides.at(i) );

                /**
                 *  Ensure that the sequence contains valid data
                 */
                char c = nucleotides.at(i);
//                if (c == '\r' || c == '\n') {
//                    continue;
//                }
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'X') {
                    cerr << "ERROR: reading sequence: " << sequence_information << endl;
                    cerr << "nucleotides: " << nucleotides << endl;
                    cerr << "nucleotide[" << i << "]: " << c << " is not 'A', 'C', 'G', 'T' or 'X'." << endl;
                    exit(0);
                }
            }
            sequences->push_back( new Sequence(sequence_information, nucleotides, max_sites, number_motifs, max_shift_distance) );
//            cout << "pushed " << "nucleotides: [" << nucleotides << "], sequence_information: [" << sequence_information << "]" << endl;

            getline(sequence_file, current_line);
            remove_carriage_return(current_line);
//            cout << "next sequence line is: " << current_line << endl;
        }

        sequence_file.close();
    } else {
        cerr << "ERROR: could not read sequence file: " << sequence_filename << endl;
        exit(0);
    }
}

void calculate_background_nucleotide_probabilities(vector<Sequence*> *sequences) {
    int total, a_count, c_count, g_count, t_count;

    total = background_pseudocounts + background_pseudocounts +  background_pseudocounts + background_pseudocounts;
    a_count = background_pseudocounts;
    c_count = background_pseudocounts;
    g_count = background_pseudocounts;
    t_count = background_pseudocounts;

    for (unsigned int i = 0; i < sequences->size(); i++) {
        Sequence* current_sequence = sequences->at(i);

        for (unsigned int j = 0; j < current_sequence->nucleotides.size(); j++) {
            switch (current_sequence->nucleotides.at(j)) {
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
            }
        }
    }

    background_nucleotide_probability[A_POSITION] = ((double)(a_count + t_count) / 2.0) / (double)total;
    background_nucleotide_probability[C_POSITION] = ((double)(c_count + g_count) / 2.0) / (double)total;
    background_nucleotide_probability[G_POSITION] = ((double)(c_count + g_count) / 2.0) / (double)total;
    background_nucleotide_probability[T_POSITION] = ((double)(a_count + t_count) / 2.0) / (double)total;

    cerr << "background probabilities:" << endl;
    cerr << "\tA:" << background_nucleotide_probability[A_POSITION] << endl;
    cerr << "\tC:" << background_nucleotide_probability[C_POSITION] << endl;
    cerr << "\tG:" << background_nucleotide_probability[G_POSITION] << endl;
    cerr << "\tT:" << background_nucleotide_probability[T_POSITION] << endl;
}


/**
 *  This counts the number of possible ways to generate 0...max_sites within a sequence
 */
void calculate_site_counts(vector<Sequence*> *sequences, vector<MotifModel> &motif_models, vector<double> &blocks, int max_sites) {
    Sequence *current_sequence;
    MotifModel *current_motif_model;

//    printf("calculating counts\n");

    for (unsigned int i = 0; i < sequences->size(); i++) {
        current_sequence = sequences->at(i);

        for (int j = 0; j < max_sites; j++) {

            for (unsigned int k = 0; k < current_sequence->nucleotides.size(); k++) {
//                printf("number_motifs [%d], current_sequence->length [%d]\n", number_motifs, current_sequence->length);

                if (k == 0) {
                    current_sequence->possible_site_counts.at(j).at(0) = 0;
                } else {
                    current_sequence->possible_site_counts.at(j).at(k) = current_sequence->possible_site_counts.at(j).at(k - 1);
                }

                for (unsigned int l = 0; l < motif_models.size(); l++) {
                    current_motif_model = &(motif_models.at(l));

                    if (current_sequence->possible_end_position(*current_motif_model, k)) {
                        //j == 0 is a special case, as all counts[-1][k] are == 1 (where -1 is 0 sites)
                        if (j == 0) {
                            if ((int)k - current_motif_model->motif_width + 1 >= 0) {

                                current_sequence->possible_site_counts.at(j).at(k) += 1;
                            }
                        } else if ((int)k - current_motif_model->motif_width >= 0) {
                            current_sequence->possible_site_counts.at(j).at(k) += current_sequence->possible_site_counts.at(j - 1).at(k - current_motif_model->motif_width);
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


        current_sequence->possible_sites.at(0) = 1 * blocks.at(0); //theres one way to have 0 sites sampled -- but we multiply this by blocks
        current_sequence->total_possible_sites = 0;

//        printf("sequence [%s]\n", current_sequence->nucleotides);
//        printf("possible_sites[0]: 1.0\n");

        for (int j = 0; j < max_sites; j++) {
            current_sequence->possible_sites.at(j+1) = blocks.at(j+1) * current_sequence->possible_site_counts.at(j).at(current_sequence->nucleotides.size() - 1);
        }

        for (int j = 0; j < max_sites + 1; j++) {
            current_sequence->total_possible_sites += current_sequence->possible_sites.at(j);
//            printf("possible_sites[%d]: %Lf\n", j, current_sequence->possible_sites[j]);
        }
//        printf("total_possible_sites: %Lf\n", current_sequence->total_possible_sites);
    }
}

long double Sequence::background_probability(double *bg_nucleotide_prob, MotifModel &motif_model, int end_position) {
    int i;
    int position;
    int nucleotide_position;
    long double probability;

    probability = 1.0;

    end_position++;

    for (i = 0; i < motif_model.motif_width; i++) {
        position = end_position - motif_model.motif_width + i;
        nucleotide_position = get_nucleotide_position( nucleotides[position] );

        probability *= bg_nucleotide_prob[nucleotide_position];
    }

    return probability;
}

long double Sequence::foreground_probability_phylogeny(MotifModel &motif_model, int end_position, PhylogenyTree *phylogeny_tree) {
    int position;
    int nucleotide_position;
    long double probability, reverse_probability;

    probability = 1.0;
    reverse_probability = 1.0;

    end_position++;
    for (int i = 0; i < motif_model.motif_width; i++) {
        position = end_position - motif_model.motif_width + i;
        nucleotide_position = get_nucleotide_position( nucleotides.at(position) );

        vector<double> *background_equilibrium = new vector<double>(background_nucleotide_probability, background_nucleotide_probability + 4);
        vector<double> *foreground_equilibrium = new vector<double>(motif_model.nucleotide_probabilities.at(i).begin(), motif_model.nucleotide_probabilities.at(i).end());
        /**
         *  
         *  Using phylogeny tree:
         *      for foreground:
         *          SubstitutionMatrix.setEquilibriumAndTTFactor( equilibrium = background [A T C G], TTFactor = 3 (maybe command line) )
         *          SubstitutionMatrix.setHB98Foreground( equilibrium_fg = [motif_model.probabilities.at(A_POSITION) motif_model.probabilities.at(T_POSITION) motif_model.probabilities.at(C_POSITION) motif_model.probabilities.at(G_POSITION)] )
         *
         *      then recursively:
         *          getSubstitutionMatrix(time_child1, response_child1);
         *          getSubstitutionMatrix(time_child2, response_child2);
         *          getSubstitutionMatrix(time_child3, response_child3);
         *          ...
         *
         *          prob_vector[i] =   (response_child1[i] * prob_vector_child1) |this is matrix multiplication| 
         *                           * (response_child2[i] * prob_vector_child2) |this is matrix multiplication|
         *                           * (response_child3[i] * prob_vector_child3) |this is matrix multiplication|
         *                           * ...
         *
         *      then:
         *          fg_prob = root_prob_vector[0] * equilibrium_fg[0] + root_prob_vector[1] * equilibrium_fg[1] + ...
         *
         *********************************
         *
         *      for bacgkround:
         *          SubstitutionMatrix.setEquilibriumAndTTFactor( equilibrium = background [A T C G], TTFactor = 3 (maybe command line) )
         *          ONLY
         *
         *      same recursion:
         *
         *      then (dot product):
         *          bg_prob = root_prob_vector[0] * background[0] + ...
         *
         */

        //must get bg_prob before fg_prob to initialize the substitution matrix correctly

        long double bg_prob = phylogeny_tree->get_background_probability(background_equilibrium, position);
        long double fg_prob = phylogeny_tree->get_foreground_probability(foreground_equilibrium);
        probability *=  fg_prob / bg_prob;

        delete background_equilibrium;
        delete foreground_equilibrium;
    }

    if (motif_model.type == MODEL_TYPE_FORWARD || motif_model.type == MODEL_TYPE_REVERSE) probability /= 2.0;

    return probability;
}

long double Sequence::foreground_probability(MotifModel &motif_model, int end_position) {
    int i;
    int position;
    int nucleotide_position;
    long double probability, reverse_probability;

    probability = 1.0;
    reverse_probability = 1.0;

    end_position++;
    for (i = 0; i < motif_model.motif_width; i++) {
        position = end_position - motif_model.motif_width + i;
        nucleotide_position = get_nucleotide_position( nucleotides.at(position) );

        probability *= motif_model.nucleotide_probabilities.at(i).at(nucleotide_position);
    }

    if (motif_model.type == MODEL_TYPE_FORWARD || motif_model.type == MODEL_TYPE_REVERSE) probability /= 2.0;

    return probability;
}

void Sequence::calculate_background_site_probabilities(vector<MotifModel> &motif_models) {
    MotifModel *current_motif;

    /**
     * Currently there's no reason we need to recalculate the background_probabilities -- this would only be done if the motifs are resized
     */
    for (unsigned int i = 0; i < motif_models.size(); i++) {
        current_motif = &(motif_models.at(i));

//        print_motif_model(current_motif);

        for (unsigned int j = 0; j < nucleotides.size(); j++) {
            if (possible_end_position(*current_motif, j)) {
                background_site_probability.at(i).at(j) = background_probability(background_nucleotide_probability, *current_motif, j);
            } else {
                background_site_probability.at(i).at(j) = 0;
            }
//          printf("%.20Lf ", sequence->background_site_probability[i][j]);
        }
    }
}

void Sequence::calculate_site_probabilities(vector<MotifModel> &motif_models, int max_sites, PhylogenyTree *phylogeny_tree) {
    long double motif_contribution;
    MotifModel *current_motif;

//    printf("calculated background probabilities\n");

    /**
     *  We always need to recalculate the site probability ratios (forground / background)
     */
//    printf("caculating site probability ratios: \n");
    for (unsigned int i = 0; i < motif_models.size(); i++) {
        current_motif = &(motif_models.at(i));

        for (unsigned int j = 0; j < nucleotides.size(); j++) {
            if (possible_end_position(*current_motif, j)) {
                if (phylogeny_tree == NULL) {
                    site_probability_ratio.at(i).at(j) = foreground_probability(*current_motif, j) / background_site_probability.at(i).at(j);
                } else {
                    site_probability_ratio.at(i).at(j) = foreground_probability_phylogeny(*current_motif, j, phylogeny_tree);
                }

                if (site_probability_ratio.at(i).at(j) < 0) {
                    cerr << endl;
                    cerr << "ERROR: calculated site probability ratio < 0: file [" << __FILE__ << "] line [" << __LINE__ << "]" << endl;
                    cout << "motif_model: " << i << endl;
                    cerr << "sequence->site_probability_ratio[" << i << "][" << j << "]: " << site_probability_ratio.at(i).at(j) << endl;
                    if (phylogeny_tree == NULL) {
                        cerr << "foreground_prob: " << foreground_probability(*current_motif, j) << ", background_prob: " << background_site_probability.at(i).at(j) << endl;
                    } else {
                        cerr << "probability_phylogeny: " << foreground_probability_phylogeny(*current_motif, j, phylogeny_tree) << endl;
                    }
                    current_motif->print(cerr);
                    exit(0);
                }
            } else {
                site_probability_ratio.at(i).at(j) = 0;
            }
//            printf("%.20Lf ", sequence->site_probability_ratio[i][j]);
        }
//        printf("\n");
    }
//    printf("calculated site probability ratios\n");

    /**
     *  Calculate the site probabilities used for sampling.
     */
    for (int j = 0; j < max_sites; j++) {
        site_probability.at(j).at(0) = 0;

        for (unsigned int k = 1; k < nucleotides.size(); k++) {
            site_probability.at(j).at(k) = site_probability.at(j).at(k - 1); //this is the probability for the background of site j,k

            for (unsigned int l = 0; l < motif_models.size(); l++) {
                current_motif = &(motif_models.at(l));

                motif_contribution = 0;
                if (possible_end_position(*current_motif, k)) {
                    if (j == 0) { 
                        //previous probability is always 1 when j == 0
                        motif_contribution = site_probability_ratio.at(l).at(k);

                        motif_probability_contribution.at(j).at(k).at(l) = motif_contribution;
                        site_probability.at(j).at(k) += motif_contribution;
                    } else if ((int)k - current_motif->motif_width < 0) {
                        //previous probability is 0
                        continue; 
                    } else {
                        motif_contribution = site_probability.at(j - 1).at(k - current_motif->motif_width) * site_probability_ratio.at(l).at(k); 
                        motif_probability_contribution.at(j).at(k).at(l) = motif_contribution;
                        site_probability.at(j).at(k) += motif_contribution;
                    }
                }

                if (site_probability.at(j).at(k) < 0) {
                    cerr << "ERROR: calculating site probabilities, probability became negative. file [" << __FILE__ << "] line [" << __LINE__ << "]" << endl;
                    cerr << "sequence->site_probability[" << j << "][" << k << "]: " << site_probability.at(j).at(k) << ", sequence->motif_probability_contribution[" << j << "][" << k << "][" << l << "]: " << motif_probability_contribution.at(j).at(k).at(l)<< endl;
                    cerr << "motif contribution: " << motif_contribution << endl;
                    if (j == 0) {
                        cerr << "j == 0" << endl;
                    } else {
                        cerr << "j == " << j << ", sequence->site_probability[" << j - 1 << "][" << k - current_motif->motif_width << "]: " << site_probability.at(j-1).at(k-current_motif->motif_width) << ", sequence->site_probability_ratio[" << l << "][" << k << "]: " << site_probability_ratio.at(l).at(k) << endl;
                        cerr << "site_probability_ratio[" << l << "][" << k-1 << "]: " << site_probability_ratio.at(l).at(k-1) << endl;
                    }
                    exit(0);
                }
            }
        }
    }
}

void Sequence::zero_accumulated_samples() {
    for (unsigned int j = 0; j < accumulated_samples.size(); j++) {
        for (unsigned int k = 0; k < accumulated_samples.at(j).size(); k++) {
            accumulated_samples.at(j).at(k) = 0;
        }
    }
}

void Sequence::print_sequence(ostream &out) {
    out << name << endl;
    out << nucleotides << endl;
}
