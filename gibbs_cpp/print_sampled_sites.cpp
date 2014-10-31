#include <iostream>
#include <fstream>
#include <vector>

#include <cstdlib>

#ifdef _LINUX64_
    #include <sys/resource.h>
#endif

#include "time.h"

#include "../../mersenne_twister/dSFMT.h"


#include "arguments.hpp"
#include "structs.hpp"
#include "sequences.hpp"
#include "motif_models.hpp"
#include "sampling.hpp"
#include "shifting.hpp"
#include "checkpoint.hpp"
#include "util.hpp"

using namespace std;

int main(int argc, char** argv) {
    int retval = 0;

    vector<string> arguments(argv, argv + argc);

    /**
     *  The following parses the values from the arguments vector.  The second argument specifies if the argument is required.
     *  If this is true and the argument is not found the program will exist with an error message. 
     */
    vector<string> motif_info;
    get_argument_vector<string>(arguments, "--motifs", true, motif_info);
    vector<MotifModel> motif_models;
    initialize_motif_models(motif_models, motif_info);

    int sample_period;
    get_argument(arguments, "--samples_period", true, sample_period);

    double best_site_percentage;
    get_argument(arguments, "--best_site_percentage", true, best_site_percentage);

    int max_sites;
    get_argument(arguments, "--max_sites", true, max_sites);

    int max_shift_distance, shift_period;
    get_arguments(arguments, "--enable_shifting", true, max_shift_distance, shift_period);

    string sequence_file;
    get_argument(arguments, "--sequence_file", true, sequence_file);

    vector<Sequence*> *sequences = new vector<Sequence*>();
    read_sequences(sequences, sequence_file, max_sites, motif_info.size() /*motif_info.size() is the number of motifs*/, max_shift_distance);

    string samples_file;
    get_argument(arguments, "--samples_file", true, samples_file);

    read_accumulated_samples(samples_file, sequences);

    double count_percentage;
    for (unsigned int j = 0; j < motif_models.size(); j++) {
        motif_models.at(j).print_short(cout);
        cout << endl;
        for (unsigned int i = 0; i < sequences->size(); i++) {

            for (unsigned int k = 0; k < sequences->at(i)->nucleotides.size(); k++) {
                count_percentage = ((double)sequences->at(i)->accumulated_samples.at(j).at(k)) / ((double)sample_period);

                if (count_percentage > best_site_percentage) {
                    printf("%5d,%2d%8d ", i, j, (k + 1) - motif_models.at(j).motif_width);
                    print_sample_and_nearest(sequences->at(i), motif_models.at(j), k, 5);
                    printf(" %7d %2.4lf %s\n", k, count_percentage, sequences->at(i)->name.c_str());
                }
            }
        }
        printf("\n");
    }

    return 0;
}

