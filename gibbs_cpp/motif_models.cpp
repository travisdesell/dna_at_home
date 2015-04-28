#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../../mersenne_twister/dSFMT.h"

#include "structs.hpp"
#include "motif_models.hpp"
#include "sequences.hpp"
#include "sampling.hpp"

using namespace std;

extern double foreground_pseudocounts;

void initialize_motif_models(vector<MotifModel> &motif_models, vector<string> &motif_info) {
    for (unsigned int i = 0; i < motif_info.size(); i++) {
        motif_models.push_back(MotifModel(motif_info.at(i)));
    }
}

MotifModel::MotifModel(string motif_info) {
    int comma_pos = motif_info.find(',');
    string type_string = motif_info.substr(0, comma_pos);

    cerr << "creating model from '" << motif_info << "'" << endl;

    string tmp = motif_info.substr(comma_pos + 1, motif_info.size() - (comma_pos + 1));

    if ( ! (stringstream(tmp) >> motif_width) ) {
        cerr << "ERROR: parsing motif model argument: '" << motif_info << endl;
        cerr << "\tcould not parse motif width." << endl;
        exit(1);
    }

    cerr << "type_string: " << type_string << ", motif_width: " << motif_width << endl;

    if (type_string.compare("normal") == 0) {
        type = MODEL_TYPE_NORMAL;
    } else if (type_string.compare("palindromic") == 0) {
        type = MODEL_TYPE_PALINDROMIC;
    } else if (type_string.compare("forward") == 0) {
        type = MODEL_TYPE_FORWARD;
    } else if (type_string.compare("reverse") == 0) {
        type = MODEL_TYPE_REVERSE;
    } else {
        cerr << "ERROR: parsing motif model argument: '" << motif_info << endl;
        cerr << "\tunknown model type: " << type_string << endl;
        exit(1);
    }

    counts = vector< vector<long> >(motif_width, vector<long>(ALPHABET_LENGTH, 0));
    nucleotide_probabilities = vector< vector<double> >(motif_width, vector<double>(ALPHABET_LENGTH, 0.0));
    reverse_nucleotide_probabilities = vector< vector<double> >(motif_width, vector<double>(ALPHABET_LENGTH, 0.0));
}

MotifModel::MotifModel(int _type, int _motif_width) {
    type = _type;
    motif_width = _motif_width;

    counts = vector< vector<long> >(motif_width, vector<long>(ALPHABET_LENGTH, 0));
    nucleotide_probabilities = vector< vector<double> >(motif_width, vector<double>(ALPHABET_LENGTH, 0.0));
    reverse_nucleotide_probabilities = vector< vector<double> >(motif_width, vector<double>(ALPHABET_LENGTH, 0.0));
}

MotifModel::MotifModel(int _type, vector< vector< double> > &_nucleotide_probabilities) {
    type = _type;
    motif_width = _nucleotide_probabilities.size();

    nucleotide_probabilities = vector< vector<double> >(_nucleotide_probabilities);

    counts = vector< vector<long> >(motif_width, vector<long>(ALPHABET_LENGTH, 0));
    reverse_nucleotide_probabilities = vector< vector<double> >(motif_width, vector<double>(ALPHABET_LENGTH, 0.0));
}

void MotifModel::zero_counts() {
    for (int j = 0; j < motif_width; j++) {
        counts.at(j).at(A_POSITION) = 0;
        counts.at(j).at(C_POSITION) = 0;
        counts.at(j).at(G_POSITION) = 0;
        counts.at(j).at(T_POSITION) = 0;
    }
}

void MotifModel::increment_counts_for_sample(Sequence *sequence, const Sample &sample) {
    int position, end_position;
    char nucleotide;

    end_position = sample.end_position + 1;

    if (end_position - motif_width < 0) {
        cerr << "ERROR (increment_counts): trying to read character with position < 0" << endl;
        cerr << "subsequence examined [" << sequence->nucleotides.substr(end_position - motif_width, motif_width) << "]" << endl;
        exit(1);
    }

    for (position = 0; position < motif_width; position++) {
        nucleotide = sequence->nucleotides.at(end_position - motif_width + position);

        switch (nucleotide) {
            case 'A': counts.at(position).at(A_POSITION)++; break;
            case 'C': counts.at(position).at(C_POSITION)++; break;
            case 'G': counts.at(position).at(G_POSITION)++; break;
            case 'T': counts.at(position).at(T_POSITION)++; break;

            default:
                    cerr << "ERROR (increment_model_counts): unknown character found in sequence " << sequence->name << "][" << sequence->nucleotides << "][" << sequence->nucleotides.size() << "] at position [" << position << "] nucleotide position [" << (end_position - motif_width + position) << "] -- [" << nucleotide << "] [" << (int)nucleotide << "]" << endl;
                    cerr << "sample->end_position [" << sample.end_position << "], motif_width [" << motif_width <<  "]" << endl;

                    cerr << "subsequence examined [" << sequence->nucleotides.substr(end_position - motif_width, motif_width) << "]" << endl;
                    cerr << "Error in file [" << __FILE__ << "], line [" << __LINE__ << "]" << endl;
        }
    }
}

void MotifModel::increment_model_counts(Sequence *sequence, const vector<Sample> &samples, int motif_model_number) {
    for (unsigned int i = 0; i < samples.size(); i++) {
        const Sample *sample = &(samples.at(i));
        if (sample->motif_model != motif_model_number) continue;

        increment_counts_for_sample(sequence, *sample);
    }
}

void increment_counts(vector<MotifModel> &motif_models, Sequence *sequence) {
    for (unsigned int i = 0; i < sequence->sampled_sites.size(); i++) {
        const Sample *sample = &(sequence->sampled_sites.at(i));

        if (sample->motif_model < 0 || sample->end_position < 0) continue;

        motif_models.at(sample->motif_model).increment_counts_for_sample(sequence, *sample);
    }
}

void decrement_counts(vector<MotifModel> &motif_models, Sequence *sequence) {
    char nucleotide;
    int position;
    int end_position;

    for (unsigned int i = 0; i < sequence->sampled_sites.size(); i++) {
        Sample sample = sequence->sampled_sites.at(i);

        if (sample.motif_model < 0 || sample.end_position < 0) continue;

        MotifModel *motif_model = &(motif_models.at(sample.motif_model));
        end_position = sample.end_position + 1;

        if (end_position - motif_model->motif_width < 0) {
            cerr << "ERROR (decrement_counts): trying to read character with position < 0" << endl;
            cerr << "subsequence examined [" << sequence->nucleotides.substr(end_position - motif_model->motif_width, motif_model->motif_width) << "]" << endl;
            exit(1);
        }

        for (position = 0; position < motif_model->motif_width; position++) {
            nucleotide = sequence->nucleotides.at(end_position - motif_model->motif_width + position);

    //        printf("incrementing count for position [%d] for nucleotide [%c]\n", position, nucleotide);

            switch (nucleotide) {
                case 'A': motif_model->counts.at(position).at(A_POSITION)--; break;
                case 'C': motif_model->counts.at(position).at(C_POSITION)--; break;
                case 'G': motif_model->counts.at(position).at(G_POSITION)--; break;
                case 'T': motif_model->counts.at(position).at(T_POSITION)--; break;

                default:
                    cerr << "ERROR (increment_model_counts): unknown character found in sequence " << sequence->name << "][" << sequence->nucleotides << "][" << sequence->nucleotides.size() << "] at position [" << position << "] nucleotide position [" << (end_position - motif_model->motif_width + position) << "] -- [" << nucleotide << "] [" << (int)nucleotide << "]" << endl;
                    cerr << "sample->end_position [" << sample.end_position << "], motif_width [" << motif_model->motif_width <<  "]" << endl;

                    cerr << "subsequence examined [" << sequence->nucleotides.substr(end_position - motif_model->motif_width, motif_model->motif_width) << "]" << endl;
                    cerr << "Error in file [" << __FILE__ << "], line [" << __LINE__ << "]" << endl;

//                    cerr << "this a possible end position?" << sequence->possible_end_position(*motif_model, position) << endl;

            }
        }
    }
}

void copy_motif_model(const MotifModel &source, MotifModel &destination) {
    int j, k;

    if (source.motif_width != destination.motif_width) {
        cerr << "ERROR copying motif model, source->width [" << source.motif_width << "] != destination->width [" << destination.motif_width << "]" << endl;
        cerr << "In file [" << __FILE__ << "], line [" << __LINE__ << "]" << endl;
        exit(1);
    }

    for (j = 0; j < source.motif_width; j++) {
        for (k = 0; k < ALPHABET_LENGTH; k++) {
            destination.nucleotide_probabilities.at(j).at(k) = source.nucleotide_probabilities.at(j).at(k);
            destination.counts.at(j).at(k) = source.counts.at(j).at(k);
        }
    }
}

void update_motif_models(vector<MotifModel> &motif_models) {
    for (unsigned int k = 0; k < motif_models.size(); k++) {
        switch (motif_models.at(k).type) {
            case MODEL_TYPE_FORWARD:
                update_motif_model_reverse_complement(motif_models.at(k), motif_models.at(k + 1));
                k++;
                break;
            default:
                motif_models.at(k).update_motif_model();
                break;
        }
    }
}

void update_motif_model_reverse_complement(MotifModel &forward_motif_model, MotifModel &reverse_motif_model) {
    double count_a, count_c, count_g, count_t;
    double sum;
    int position;

    switch (forward_motif_model.type) {
        case MODEL_TYPE_FORWARD:
                for (position = 0; position < forward_motif_model.motif_width; position++) {

                    int reverse_pos = forward_motif_model.motif_width - (position + 1);

                    count_a = foreground_pseudocounts + forward_motif_model.counts.at(position).at(A_POSITION) + reverse_motif_model.counts.at(reverse_pos).at(A_POSITION_REVERSE);
                    count_c = foreground_pseudocounts + forward_motif_model.counts.at(position).at(C_POSITION) + reverse_motif_model.counts.at(reverse_pos).at(C_POSITION_REVERSE);
                    count_g = foreground_pseudocounts + forward_motif_model.counts.at(position).at(G_POSITION) + reverse_motif_model.counts.at(reverse_pos).at(G_POSITION_REVERSE);
                    count_t = foreground_pseudocounts + forward_motif_model.counts.at(position).at(T_POSITION) + reverse_motif_model.counts.at(reverse_pos).at(T_POSITION_REVERSE);

                    sum = count_a + count_c + count_g + count_t;

                    forward_motif_model.nucleotide_probabilities.at(position).at(A_POSITION) = count_a / sum; 
                    forward_motif_model.nucleotide_probabilities.at(position).at(C_POSITION) = count_c / sum; 
                    forward_motif_model.nucleotide_probabilities.at(position).at(G_POSITION) = count_g / sum; 
                    forward_motif_model.nucleotide_probabilities.at(position).at(T_POSITION) = count_t / sum; 

                    reverse_motif_model.nucleotide_probabilities.at(reverse_pos).at(A_POSITION_REVERSE) = forward_motif_model.nucleotide_probabilities.at(position).at(A_POSITION);
                    reverse_motif_model.nucleotide_probabilities.at(reverse_pos).at(C_POSITION_REVERSE) = forward_motif_model.nucleotide_probabilities.at(position).at(C_POSITION);
                    reverse_motif_model.nucleotide_probabilities.at(reverse_pos).at(G_POSITION_REVERSE) = forward_motif_model.nucleotide_probabilities.at(position).at(G_POSITION);
                    reverse_motif_model.nucleotide_probabilities.at(reverse_pos).at(T_POSITION_REVERSE) = forward_motif_model.nucleotide_probabilities.at(position).at(T_POSITION);
                }
                break;

        default:
                fprintf(stderr, "ERROR: unknown model type in update_motif_model_reverse_complement, FILE [%s], LINE [%d]\n", __FILE__, __LINE__);
                fprintf(stderr, "FORWARD MODEL:\n");
                forward_motif_model.print(cerr);

                fprintf(stderr, "\n");
                fprintf(stderr, "REVERSEMODEL:\n");
                reverse_motif_model.print(cerr);
                exit(1);
    }

}

void MotifModel::update_motif_model() {
    double count_a, count_c, count_g, count_t;
    double sum;
    int position;

    switch (type) {
        case MODEL_TYPE_NORMAL:
                for (position = 0; position < motif_width; position++) {

                    count_a = foreground_pseudocounts + counts.at(position).at(A_POSITION);
                    count_c = foreground_pseudocounts + counts.at(position).at(C_POSITION);
                    count_g = foreground_pseudocounts + counts.at(position).at(G_POSITION);
                    count_t = foreground_pseudocounts + counts.at(position).at(T_POSITION);

                    sum = count_a + count_c + count_g + count_t;

                    nucleotide_probabilities.at(position).at(A_POSITION) = count_a / sum; 
                    nucleotide_probabilities.at(position).at(C_POSITION) = count_c / sum; 
                    nucleotide_probabilities.at(position).at(G_POSITION) = count_g / sum; 
                    nucleotide_probabilities.at(position).at(T_POSITION) = count_t / sum; 

            //        printf("set probs[%d] to %10.5lf %10.5lf %10.5lf %10.5lf\n", position, (count_a + 0.28) / sum, (count_c + 0.28) / sum, (count_g + 0.28) / sum, (count_t + 0.28) / sum);
                }
                break;

        case MODEL_TYPE_PALINDROMIC:
                for (position = 0; position < (motif_width / 2); position++) {
//                    printf("calculating palindromic probability for position [%d] and [%d]\n", position, (motif_width - (position + 1)));

//                    printf("foreground_pseudocounts: %lf\n", foreground_pseudocounts);
                    int pos = motif_width - (position + 1);
                    count_a = foreground_pseudocounts + counts.at(position).at(A_POSITION) + counts.at(pos).at(A_POSITION_REVERSE);
                    count_c = foreground_pseudocounts + counts.at(position).at(C_POSITION) + counts.at(pos).at(C_POSITION_REVERSE);
                    count_g = foreground_pseudocounts + counts.at(position).at(G_POSITION) + counts.at(pos).at(G_POSITION_REVERSE);
                    count_t = foreground_pseudocounts + counts.at(position).at(T_POSITION) + counts.at(pos).at(T_POSITION_REVERSE);

                    sum = count_a + count_c + count_g + count_t;

                    nucleotide_probabilities.at(position).at(A_POSITION) = count_a / sum; 
                    nucleotide_probabilities.at(position).at(C_POSITION) = count_c / sum; 
                    nucleotide_probabilities.at(position).at(G_POSITION) = count_g / sum; 
                    nucleotide_probabilities.at(position).at(T_POSITION) = count_t / sum; 

                    nucleotide_probabilities.at(pos).at(A_POSITION_REVERSE) = count_a / sum; 
                    nucleotide_probabilities.at(pos).at(C_POSITION_REVERSE) = count_c / sum; 
                    nucleotide_probabilities.at(pos).at(G_POSITION_REVERSE) = count_g / sum; 
                    nucleotide_probabilities.at(pos).at(T_POSITION_REVERSE) = count_t / sum; 

////                    printf("probabilities: [%lf][%lf][%lf][%lf]\n", motif_model->nucleotide_probabilities[position][A_POSITION], motif_model->nucleotide_probabilities[position][C_POSITION], motif_model->nucleotide_probabilities[position][G_POSITION], motif_model->nucleotide_probabilities[position][T_POSITION]);

//                    printf("probabilities: [%lf][%lf][%lf][%lf]\n", motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][A_POSITION], motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][C_POSITION], motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][G_POSITION], motif_model->nucleotide_probabilities[motif_model->motif_width - (position + 1)][T_POSITION]);
                }

                if ((motif_width % 2) == 1) {
                    position = (motif_width / 2) + 1;
//                    printf("palindromic middle for model width [%d] is [%d]\n", motif_model->motif_width, position);

                    count_a = foreground_pseudocounts + counts.at(position).at(A_POSITION);
                    count_c = foreground_pseudocounts + counts.at(position).at(C_POSITION);
                    count_g = foreground_pseudocounts + counts.at(position).at(G_POSITION);
                    count_t = foreground_pseudocounts + counts.at(position).at(T_POSITION);

                    sum = count_a + count_c + count_g + count_t;

                    nucleotide_probabilities.at(position).at(A_POSITION) = count_a / sum; 
                    nucleotide_probabilities.at(position).at(C_POSITION) = count_c / sum; 
                    nucleotide_probabilities.at(position).at(G_POSITION) = count_g / sum; 
                    nucleotide_probabilities.at(position).at(T_POSITION) = count_t / sum; 
                }

                break;
        default:
                cerr << "ERROR: unknown model type in update_motif_model, FILE [" << __FILE__ << "], LINE [" << __LINE__ << "]" << endl;
                print(cerr);
                exit(1);
    }
}

void MotifModel::print_short(ostream &out_stream) {
    out_stream << MODEL_TYPES[type] << "," << motif_width;
}

void MotifModel::print(ostream &out_stream) {
    int j, k;

    out_stream << "Model type: " << MODEL_TYPES[type] << endl;
    out_stream << "Model length: " << motif_width << endl;

    for (j = 0; j < ALPHABET_LENGTH; j++) {
        out_stream << "letter count [" << ALPHABET[j] << "]: ";

        for (k = 0; k < motif_width; k++) {
            out_stream << " " << counts.at(k).at(j);
        }
        out_stream << endl;
    }
    out_stream << endl;

    for (j = 0; j < ALPHABET_LENGTH; j++) {
        out_stream << "letter prob [" << ALPHABET[j] << "]: ";

        for (k = 0; k < motif_width; k++) {
            out_stream << " " << nucleotide_probabilities.at(k).at(j);
        }
        out_stream << endl;
    }
    out_stream << endl;
}

string MotifModel::generate_possible_nucleotides() {
    string s = "";

    for (int i = 0; i < motif_width; i++) {
        double value = dsfmt_gv_genrand_open_close();

        if (value < nucleotide_probabilities.at(i).at(A_POSITION)) {
            s.append("A");
        } else {
            value -= nucleotide_probabilities.at(i).at(A_POSITION);
            if (value < nucleotide_probabilities.at(i).at(C_POSITION)) {
                s.append("C");
            } else {
                value -= nucleotide_probabilities.at(i).at(C_POSITION);
                if (value < nucleotide_probabilities.at(i).at(G_POSITION)) {
                    s.append("G");
                } else {
                    value -= nucleotide_probabilities.at(i).at(G_POSITION);
                    if (value < nucleotide_probabilities.at(i).at(T_POSITION)) {
                        s.append("T");
                    } else {
                        cerr << "ERROR(" << __FILE__ << "," << __LINE__ << "): could not generate possible nucleotides for motif: " << endl;
                        print(cerr);
                        cerr << endl << "random value to choose nucleotide " << i << " did not select any random nucleotide (chances are the nucleotide probabilities did not sum to 1" << endl;
                        exit(1);
                    }                  
                }                  
            } 
        }
    }

    return s;
}

void read_motifs_from_file(string filename, vector<MotifModel> &motif_models) {
    ifstream in(filename.c_str());

    int type = -1;

    while (in.good() && !in.eof()) {
        bool is_reverse_complement = false;
        string line;
        getline(in, line);
//        cout << "type line: " << line << endl;

        if (line.compare(">reverse_complement") == 0) {
            is_reverse_complement = true;
            type = MODEL_TYPE_FORWARD;
        } else if (line.compare(">palindromic") == 0) {
            type = MODEL_TYPE_PALINDROMIC;
        } else if (line.compare(">forward") == 0) {
            type = MODEL_TYPE_FORWARD;
        } else if (line.compare(">reverse") == 0) {
            type = MODEL_TYPE_REVERSE;
        }

        vector< vector<double> > nucleotide_probabilities;
        vector< vector<double> > reverse_nucleotide_probabilities;

        int previous_width = 0;
        for (int i = 0; i < ALPHABET_LENGTH; i++) {
            getline(in, line);
//            cout << "read line: " << line << endl;
            istringstream iss(line);

            int current = 0;
            while (iss.good() && !iss.eof()) {
                if (i == 0) nucleotide_probabilities.push_back( vector<double>(ALPHABET_LENGTH, 0.0) );

                double value;
                iss >> value;
//                cout << "value: " << value << endl;
                
                nucleotide_probabilities.at(current).at(i) = value;
                current++;
            }
            if (i > 0) {
                if (current != previous_width) {
                    cerr << "ERROR(" << __FILE__ << ", " << __LINE__ << "): reading motif from file, different widths for nucleotides" << endl;
                    cerr << "  current line: '" << line << "' has " << current << " nucleotides." << endl;
                    cerr << "  previous line had " << previous_width << " nucleotides." << endl;
                    exit(1);
                }
            }
            previous_width = current;
        }

        motif_models.push_back( MotifModel(type, nucleotide_probabilities) );
//        motif_models.back().print(cout);

        if (is_reverse_complement) {
            reverse_nucleotide_probabilities = vector< vector<double> >(nucleotide_probabilities.size(), vector<double>(ALPHABET_LENGTH, 0.0));

            int i_size = nucleotide_probabilities.size();
            for (int i = i_size - 1; i >= 0; i--) {
                int j_size = nucleotide_probabilities.at(i).size();
                for (int j = j_size - 1; j >= 0; j--) {
                    int rev_i = (i_size - 1) - i;
                    int rev_j = (j_size - 1) - j;
//                    cout << "setting " << i << "," << j << " set to " << rev_i << "," << rev_j << endl;
                    reverse_nucleotide_probabilities.at(i).at(j) = nucleotide_probabilities.at(rev_i).at(rev_j);
                }
            }

            motif_models.push_back( MotifModel(MODEL_TYPE_REVERSE, reverse_nucleotide_probabilities) );
//            motif_models.back().print(cout);
        }

        getline(in, line);
    }
}
