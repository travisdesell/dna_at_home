#include <fstream>
using std::ifstream;
using std::ostream;

#include <iomanip>
using std::fixed;
using std::setw;
using std::setprecision;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <sstream>
using std::stringstream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "mpi.h"

//from TAO
#include "mpi/mpi_particle_swarm.hxx"
#include "mpi/mpi_differential_evolution.hxx"

#include "asynchronous_algorithms/particle_swarm.hxx"
#include "asynchronous_algorithms/differential_evolution.hxx"

//from undvc_common
#include "arguments.hxx"

void split_line_by(string line, char seperator, vector<string> &result) {
    stringstream test(line);
    string segment;

    //parse the line by commas
    int count = 0;
    while(getline(test, segment, seperator)) {
        //cout << "pushing to row[" << count << "]: '" << segment << "'" << endl;
        count++;

        result.push_back(segment);
    }   
}

int is_valid_nucleotide(char nucleotide) {
    switch (nucleotide) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'X':
            return 1;
        default:
            return 0;
    }
}

class PositionWeightMatrix {
    private:
        double fitness;

        //index by A C G T
        vector< vector<double> > nucleotide_probabilities;

    public:
        double get_fitness() const {
            return fitness;
        }

        void set_fitness(double _fitness) {
            fitness = _fitness;
        }

        int model_length() const {
            return nucleotide_probabilities.size();
        }

        double probability(int position, int nucleotide) const {
            return nucleotide_probabilities.at(position).at(nucleotide);
        }

        PositionWeightMatrix() {
            fitness = 0.0;
        }

        PositionWeightMatrix(const vector<double> &counts) {
            fitness = 0.0;
            nucleotide_probabilities = vector< vector<double> >(counts.size() / 4, vector<double>(4, 0.0));

            for (uint32_t i = 0; i < nucleotide_probabilities.size(); i++) {
                nucleotide_probabilities[i][0] = counts[i * 4];
                nucleotide_probabilities[i][1] = counts[(i * 4) + 1];
                nucleotide_probabilities[i][2] = counts[(i * 4) + 2];
                nucleotide_probabilities[i][3] = counts[(i * 4) + 3];

                double sum = nucleotide_probabilities[i][0] + nucleotide_probabilities[i][1] + nucleotide_probabilities[i][2] + nucleotide_probabilities[i][3];

                nucleotide_probabilities[i][0] /= sum;
                nucleotide_probabilities[i][1] /= sum;
                nucleotide_probabilities[i][2] /= sum;
                nucleotide_probabilities[i][3] /= sum;
            }
        }

        PositionWeightMatrix(const vector< vector<double> > &_nucleotide_probabilities) {
            fitness = 0.0;
            nucleotide_probabilities = _nucleotide_probabilities;
        }

        PositionWeightMatrix(string pwm_filename) {
            fitness = 0.0;

            cerr << "opening position weight matrix file '" << pwm_filename << "' for reading." << endl;
            ifstream pwm_file(pwm_filename.c_str());

            string line;

            while (getline(pwm_file, line)) {
                //erase all carraige returns from the line (thanks windows/excel)
                line.erase(remove(line.begin(), line.end(), '\r'), line.end());

                //cout << "line: '" << line << "'" << endl;

                if (line.size() == 0 || line[0] == '#') continue;

                vector<string> row;
                split_line_by(line, ',', row);

                if (row.size() != 4) {
                    cerr << "Error reading line from PWM file: '" << line << "'" << endl;
                    cerr << "Did not have 4 nucelotide weights." << endl;
                    exit(1);
                }

                vector<double> probabilities;
                probabilities.push_back(stof(row[0]));
                probabilities.push_back(stof(row[1]));
                probabilities.push_back(stof(row[2]));
                probabilities.push_back(stof(row[3]));

                double sum = probabilities[0] + probabilities[1] + probabilities[2] + probabilities[3];

                probabilities[0] /= sum;
                probabilities[1] /= sum;
                probabilities[2] /= sum;
                probabilities[3] /= sum;

                nucleotide_probabilities.push_back(probabilities);
            }

            pwm_file.close();
        }

        friend ostream &operator<<(std::ostream &os, PositionWeightMatrix const &pwm);
};

ostream &operator<<(ostream &os, PositionWeightMatrix const &pwm) { 
    os << "position weight matrix:" << endl;
    os << "\t" << setw(5) << "A";
    os << ", " << setw(5) << "C";
    os << ", " << setw(5) << "G";
    os << ", " << setw(5) << "T";
    os << endl;

    string motif_nucleotides("");
    for (uint32_t i = 0; i < pwm.model_length(); i++) {
        os << "\t" << fixed << setw(5) << setprecision(4) << pwm.probability(i, 0);
        os << ", " << fixed << setw(5) << setprecision(4) << pwm.probability(i, 1);
        os << ", " << fixed << setw(5) << setprecision(4) << pwm.probability(i, 2);
        os << ", " << fixed << setw(5) << setprecision(4) << pwm.probability(i, 3);
        os << endl;

        if (pwm.probability(i, 0) > 0.90) motif_nucleotides.append("A");
        else if (pwm.probability(i, 1) > 0.90) motif_nucleotides.append("C");
        else if (pwm.probability(i, 2) > 0.90) motif_nucleotides.append("G");
        else if (pwm.probability(i, 3) > 0.90) motif_nucleotides.append("T");
        else motif_nucleotides.append("?");
    }

    os << "motif: " << motif_nucleotides << endl;

    return os; 
}


class Sequence {
    private:
        uint32_t id;
        string name;
        string nucleotides;

    public:
        string get_name() const {
            return name;
        }

        string get_nucleotides() const {
            return nucleotides;
        }

        Sequence(uint32_t _id, string _name, string _nucleotides) {
            id = _id;
            name = _name;
            nucleotides = _nucleotides;
        }

        double calculate_fitness(const vector<PositionWeightMatrix> &pwms, bool verbose) const {
            double total_fitness = 0.0;

            double model_length = pwms[0].model_length();

            for (uint32_t i = 0; i < nucleotides.size() - model_length; i++) {
                vector<double> pwm_forward_fitness(pwms.size(), 1.0);
                vector<double> pwm_reverse_fitness(pwms.size(), 1.0);


                for (uint32_t k = 0; k < pwms.size(); k++) {
                    for (uint32_t j = 0; j < model_length; j++) {
                        switch (nucleotides[i + j]) {
                            case 'A': pwm_forward_fitness[k] *= pwms[k].probability(j, 0);
                                      break;
                            case 'C': pwm_forward_fitness[k] *= pwms[k].probability(j, 1);
                                      break;
                            case 'G': pwm_forward_fitness[k] *= pwms[k].probability(j, 2);
                                      break;
                            case 'T': pwm_forward_fitness[k] *= pwms[k].probability(j, 3);
                                      break;
                            case 'X': pwm_forward_fitness[k] = 0;
                                      break;
                        }

                        switch(nucleotides[i + (6 - j)]) {
                            case 'A': pwm_reverse_fitness[k] *= pwms[k].probability(j, 0);
                                      break;
                            case 'C': pwm_reverse_fitness[k] *= pwms[k].probability(j, 1);
                                      break;
                            case 'G': pwm_reverse_fitness[k] *= pwms[k].probability(j, 2);
                                      break;
                            case 'T': pwm_reverse_fitness[k] *= pwms[k].probability(j, 3);
                                      break;
                            case 'X': pwm_reverse_fitness[k] = 0;
                                      break;
                        }
                    }

                    if (verbose) {
                        if (pwm_forward_fitness[k] > 0.5) {
                            cout << "high fitness: " << pwm_forward_fitness[k] << " for sequence " << id << " position " << i << " '";
                            for (uint32_t j = 0; j < model_length; j++) {
                                cout << nucleotides[i + j];
                            }
                            cout << "'" << endl;
                        }

                        if (pwm_reverse_fitness[k] > 0.5) {
                            cout << "high fitness: " << pwm_reverse_fitness[k] << " for sequence " << id << " reverse complement position " << i << " '";
                            for (uint32_t j = 0; j < model_length; j++) {
                                switch (nucleotides[i + (6 - j)]) {
                                    case 'A': cout << "T"; break;
                                    case 'C': cout << "G"; break;
                                    case 'G': cout << "C"; break;
                                    case 'T': cout << "A"; break;
                                    case 'X': cout << "X"; break;
                                }
                            }
                            cout << "'" << endl;
                        }
                    }
                }

                double fitness_value = 1.0;
                double final_fitness_value = 0.0;
                for (uint32_t k = 0; k < pwms.size(); k++) {
                    if (pwm_forward_fitness[k] > pwm_reverse_fitness[k]) {
                        pwm_forward_fitness[k] = fitness_value * pwm_forward_fitness[k];
                        fitness_value = fitness_value - pwm_forward_fitness[k];
                    } else {
                        pwm_reverse_fitness[k] = fitness_value * pwm_reverse_fitness[k];
                        fitness_value = fitness_value - pwm_reverse_fitness[k];
                    }

                    double max_A_fitness = 0.0;
                    double max_C_fitness = 0.0;
                    double max_G_fitness = 0.0;
                    double max_T_fitness = 0.0;

                    for (uint32_t j = 0; j < model_length; j++) {
                        if (max_A_fitness < pwms[k].probability(j, 0)) max_A_fitness = pwms[k].probability(j,0);
                        if (max_C_fitness < pwms[k].probability(j, 1)) max_C_fitness = pwms[k].probability(j,1);
                        if (max_G_fitness < pwms[k].probability(j, 2)) max_G_fitness = pwms[k].probability(j,2);
                        if (max_T_fitness < pwms[k].probability(j, 3)) max_T_fitness = pwms[k].probability(j,3);
                    }

                    double fitness_mod = max_A_fitness * max_C_fitness * max_G_fitness * max_T_fitness;
                    //double fitness_mod = (max_A_fitness + max_C_fitness + max_G_fitness + max_T_fitness + 0.4) / 4.4;

                    pwm_forward_fitness[k] *= fitness_mod;
                    pwm_reverse_fitness[k] *= fitness_mod;

                    if (pwm_forward_fitness[k] > pwm_reverse_fitness[k]) {
                        final_fitness_value += pwm_forward_fitness[k];
                    } else {
                        final_fitness_value += pwm_reverse_fitness[k];
                    }
                }

                /*
                if (fitness_value < 0.01) {
                    fitness_value = 1.0;
                    for (uint32_t k = 0; k < pwms.size(); k++) {
                        if (pwm_forward_fitness[k] > pwm_reverse_fitness[k]) {
                            cout << "fitness_value = " << fitness_value << " - (" << fitness_value << " * " << pwm_forward_fitness[k] << ") -- forward, reverse was: " << pwm_reverse_fitness[k] << endl;
                            fitness_value = fitness_value - (fitness_value * pwm_forward_fitness[k]);
                        } else {
                            cout << "fitness_value = " << fitness_value << " - (" << fitness_value << " * " << pwm_reverse_fitness[k] << ") -- reverse, forward was: " << pwm_forward_fitness[k] << endl;
                            fitness_value = fitness_value - (fitness_value * pwm_reverse_fitness[k]);
                        }
                    }

                    cout << "final fitness: " << (1.0 - fitness_value) << endl << endl;
                }
                */

                //total_fitness += 1.0 - fitness_value;
                total_fitness += final_fitness_value;
            }

            return total_fitness;
        }


        friend ostream &operator<<(std::ostream &os, Sequence const &sequence);
};

ostream &operator<<(ostream &os, Sequence const &sequence) { 
    os << sequence.get_name() << endl;
    os << sequence.get_nucleotides() << endl;
    return os; 
}

class Sequences {
    private:
        bool verbose;
        vector<Sequence> sequences;

    public:

        uint32_t number_sequences() const {
            return sequences.size();
        }

        Sequence get(uint32_t i) const {
            return sequences.at(i);
        }

        Sequences() {
            verbose = true;
        }

        void set_verbose(bool _verbose) {
            verbose = _verbose;
        }

        void read_from_file(string sequences_filename) {
            cerr << "opening sequences file '" << sequences_filename << "' for reading." << endl;
            ifstream sequence_file(sequences_filename.c_str());

            string line;

            uint32_t id = 0;
            while (getline(sequence_file, line)) {
                //erase all carraige returns from the line (thanks windows/excel)
                line.erase(remove(line.begin(), line.end(), '\r'), line.end());


                /**
                 *  File format is (for each sequence):
                 *  >sequence information
                 *  NUCLEOTIDES
                 *  NUCLEOTIDES
                 *  NUCLEOTIDES
                 *  <blank line>
                 */
                string sequence_information = line;
                //cout << "sequence_information: '" << sequence_information << "'" << endl;

                string nucleotides;
                getline(sequence_file, line);
                line.erase(remove(line.begin(), line.end(), '\r'), line.end());
                //cout << "current line: '" << line << "'" << endl;

                do {
                    nucleotides.append(line);
                    //cout << "nucleotides: " << nucleotides << endl;

                    getline(sequence_file, line);
                    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
                    //cout << "current line: '" << line << "'" << endl;
                } while (sequence_file.good() && !sequence_file.eof() && line.compare("") != 0);

                for (unsigned int i = 0; i < nucleotides.size(); i++) {
                    /**
                     *  Convert the sequence to all uppercase
                     */
                    nucleotides.at(i) = toupper( (unsigned char)nucleotides.at(i) );

                    /**
                     *  Ensure that the sequence contains valid data
                     */
                    if (!is_valid_nucleotide(nucleotides.at(i))) {
                        cerr << "ERROR: reading sequence: " << sequence_information << endl;
                        cerr << "nucleotides: " << nucleotides << endl;
                        cerr << "nucleotide[" << i << "]: " << nucleotides.at(i) << " is not 'A', 'C', 'G', 'T' or 'X'." << endl;
                        exit(1);
                    }
                }
                sequences.push_back(Sequence(id, sequence_information, nucleotides));
                id++;

                //cout << "pushed back sequence, getting next line!" << endl << endl;
            }

            sequence_file.close();
        }

        /*
        double calculate_fitness(const PositionWeightMatrix &pwm) const {
            double fitness = 0.0;

            for (uint32_t i = 0; i < sequences.size(); i++) {
                fitness += sequences.at(i).calculate_fitness(pwm, verbose);
            }

            return fitness;
        }
        */

        double calculate_fitness(const vector<PositionWeightMatrix> &pwms) const {
            double fitness = 0.0;

            for (uint32_t i = 0; i < sequences.size(); i++) {
                fitness += sequences.at(i).calculate_fitness(pwms, verbose);
            }

            return fitness;
        }


        friend ostream &operator<<(std::ostream &os, Sequences const &sequences);
};

ostream &operator<<(ostream &os, Sequences const &sequences) { 
    os << "SEQUENCES:" << endl;
    for (uint32_t i = 0; i < sequences.number_sequences(); i++) {
        os << sequences.get(i) << endl;
    }
    os << endl;
    return os; 
}

uint32_t number_motifs;
uint32_t parameters_per_motif;
Sequences sequences;

double objective_function(const vector<double> &parameters) {
    vector<PositionWeightMatrix> pwms(number_motifs);
    for (uint32_t i = 0; i < number_motifs; i++) {
        auto start = parameters.begin() + (i *  parameters_per_motif);
        auto end = start + parameters_per_motif;
        vector<double> motif_parameters(start, end);

        pwms[i] = PositionWeightMatrix(motif_parameters);
    }

    //cout << "evaluating " << pvm << endl;
    //cout << sequences << endl;

    double fitness = sequences.calculate_fitness(pwms);

    //cout << "calculated fitness: " << fitness << endl;

    return fitness;
}

struct sort_pwms_desc {
    bool operator()(const PositionWeightMatrix &pwm1, const PositionWeightMatrix &pwm2) {
        return pwm1.get_fitness() > pwm2.get_fitness();
    }   
};


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, max_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &max_rank);

    vector<string> arguments(argv, argv + argc);

    string sequences_filename = arguments[1];
    string position_weight_matrix_filename = arguments[2];

    sequences.read_from_file(sequences_filename);
    sequences.set_verbose(false);

    //cout << sequences << endl;

    PositionWeightMatrix pwm(position_weight_matrix_filename);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        cout << pwm << endl;

        vector<PositionWeightMatrix> pwms;
        pwms.push_back(pwm);
        cout << "fitness: " << sequences.calculate_fitness(pwms) << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    get_argument(arguments, "--number_motifs", true, number_motifs);
    parameters_per_motif = pwm.model_length() * 4;

    vector<double> min_bound(pwm.model_length() * 4 * number_motifs, 0.0);
    vector<double> max_bound(pwm.model_length() * 4 * number_motifs, 5.0);

    vector<double> best;

    string search_type;
    get_argument(arguments, "--search_type", true, search_type);
    if (search_type.compare("ps") == 0) {
        ParticleSwarm ps(min_bound, max_bound, arguments);
        ps.iterate(objective_function);
        best = ps.get_global_best();

    } else if (search_type.compare("de") == 0) {
        DifferentialEvolution de(min_bound, max_bound, arguments);
        de.iterate(objective_function);
        best = de.get_global_best();

    } else if (search_type.compare("de_mpi") == 0) {
        DifferentialEvolutionMPI de(min_bound, max_bound, arguments);
        de.go(objective_function);
        best = de.get_global_best();

    } else if (search_type.compare("ps_mpi") == 0) {
        ParticleSwarmMPI ps(min_bound, max_bound, arguments);
        ps.go(objective_function);
        best = ps.get_global_best();


    } else {
        cerr << "Improperly specified search type: '" << search_type.c_str() <<"'" << endl;
        cerr << "Possibilities are:" << endl;
        cerr << "    de     -       differential evolution" << endl;
        cerr << "    ps     -       particle swarm optimization" << endl;
        exit(1);

    }

    //sequences.set_verbose(true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        cout << "Best PWMs:" << endl;

        vector<PositionWeightMatrix> all_pwms;

        for (uint32_t i = 0; i < number_motifs; i++) {
            vector<PositionWeightMatrix> pwms(1);

            auto start = best.begin() + (i *  parameters_per_motif);
            auto end = start + parameters_per_motif;
            vector<double> motif_parameters(start, end);

            pwms[0] = PositionWeightMatrix(motif_parameters);
            pwms[0].set_fitness( sequences.calculate_fitness(pwms) );

            all_pwms.push_back(pwms[0]);

            //cout << pwms[0] << endl;
            //cout << "fitness: " << sequences.calculate_fitness(pwms) << endl;
        }

        sort(all_pwms.begin(), all_pwms.end(), sort_pwms_desc());

        cout << endl;
        cout << "SORTED BEST PWMS:" << endl;
        for (uint32_t i = 0; i < all_pwms.size(); i++) {
            cout << all_pwms[i];
            cout << "fitness: " << all_pwms[i].get_fitness() << endl << endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

//    double best_fitness = objective_function(ps.get_global_best());
//    cout << "fitness: " << best_fitness << endl;
}
