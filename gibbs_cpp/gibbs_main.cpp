#include <iostream>
#include <vector>

#ifdef _WIN32
    #include "boinc_win.h"
    #include "str_util.h"
#endif

#ifdef _BOINC_
    #include "diagnostics.h"
    #include "util.h"
    #include "filesys.h"
    #include "boinc_api.h"
    #include "mfile.h"
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


#define SITES_CHECKPOINT_FILE "gibbs_sites_checkpoint.txt"
#define SAMPLES_CHECKPOINT_FILE "gibbs_samples_checkpoint.txt"
#define SAMPLES_OUTPUT_FILE "accumulated_samples.txt"

using namespace std;

int max_sites;

int number_sequences;
int max_shift_distance;

double background_nucleotide_probability[ALPHABET_LENGTH];

double foreground_pseudocounts = 0.28;
double background_pseudocounts = 5.0;


int number_motifs;
MotifModel **motif_models;

int main(int argc, char** argv) {
    int retval;
    double count_percentage;
    double progress;

    retval = 0;
    #ifdef _BOINC_
        #ifdef BOINC_APP_GRAPHICS
            #if defined(_WIN32) || defined(__APPLE)
                retval = boinc_init_graphics(worker);
            #else
                retval = boinc_init_graphics(worker, argv[0]);
            #endif
        #else
            retval = boinc_init();
        #endif

        if (retval) exit(retval);
    #endif

    int independent_walk = 0;
    int total_independent_walks = 1;
    int seed = -1;
    int max_shift_distance = -1;
    int shift_period = -1;
    double print_best_sites = 0.0;
    int print_current_sites_frequency = 0;
    int max_sites = -1;
    int iteration = 0;
    int sample_period = 0;
    int burn_in_period = 0;
    int sites_from_arguments = 0;

    vector<string> arguments(argv, argv + argc);

    string sequence_file;
    string sites_file;
    vector<string> motif_info;

    /**
     *  The following parses the values from the arguments vector.  The second argument specifies if the argument is required.
     *  If this is true and the argument is not found the program will exist with an error message. 
     */
    get_argument_vector<string>(arguments, "--motifs", true, motif_info);
    get_arguments(arguments, "--enable_shifting", false, max_shift_distance, shift_period);

    get_argument(arguments, "--sequence_file", true, sequence_file);
    get_argument(arguments, "--current_sites", false, sites_file);
    get_argument(arguments, "--burn_in_period", true, burn_in_period);
    get_argument(arguments, "--sample_period", false, sample_period);
    get_argument(arguments, "--max_sites", true, max_sites);
    get_argument(arguments, "--seed", false, seed);
    get_argument(arguments, "--print_best_sites", false, print_best_sites);
    get_argument(arguments, "--total_independent_walks", false, total_independent_walks);
    get_argument(arguments, "--print_current_sites_frequency", false, print_current_sites_frequency);
    int print_current_sites_frequency_initial = print_current_sites_frequency;

    bool print_current_sites                = argument_exists(arguments, "--print_current_sites");
    bool print_accumulated_samples          = argument_exists(arguments, "--print_accumulated_samples");
    bool print_motif_models                 = argument_exists(arguments, "--print_motif_models");
    bool print_current_sites_logarithmic    = argument_exists(arguments, "--print_current_sites_logarithmic");
    bool use_checkpointing                  = !argument_exists(arguments, "--no_checkpointing");

    vector<double> blocks;
    if (argument_exists(arguments, "--blocks")) {
        get_argument_vector<double>(arguments, "--blocks", false, blocks);

        if ((int)blocks.size() != max_sites + 1) {
            cerr << "ERROR: number of blocks (--blocks <b1> <b2> ... <bn>) not equal to max_sites + 1 (--max_sites <n>)" << endl;
            cerr << "\tblocks:";
            for (unsigned int i = 0; i < blocks.size(); i++) cerr << " " << blocks.at(i);
            cerr << endl;
            cerr << "\tmax_sites + 1: " << max_sites + 1 << endl;
            exit(0);
        }
    } else {
        for (int j = 0; j < max_sites + 1; j++) blocks.push_back(1.0 / ((double)max_sites + 1.0));
    }

    cerr << "blocks:";
    for (unsigned int i = 0; i < blocks.size(); i++) cerr << " " << blocks.at(i);
    cerr << endl;

    vector<Sequence> sequences;
    read_sequences(sequences, sequence_file, max_sites, motif_info.size() /*motif_info.size() is the number of motifs*/);
    /**
     *  Calculate the background probabilities since we only need to do this once.
     */
    calculate_background_nucleotide_probabilities(sequences);
    if (sites_file.compare("") != 0) read_sites(sites_file, sequences);

    vector<MotifModel> motif_models;
    initialize_motif_models(motif_models, motif_info);

    int starting_from_checkpoint = 0;
    if (use_checkpointing) {
        string cp_file = SITES_CHECKPOINT_FILE;
        starting_from_checkpoint = read_sites_from_checkpoint(cp_file, sequences, seed, iteration, independent_walk);
        if (starting_from_checkpoint == 0) {
            cerr << "not starting from checkpoint" << endl;
            //Sample from the sequences uniformly at random if there is no checkpoint file.
        } else {
            cerr << "starting from checkpoint" << endl;
            starting_from_checkpoint = 1;
            if (iteration >= burn_in_period) {
                read_accumulated_samples(string(SAMPLES_CHECKPOINT_FILE), sequences);
            }
        }
    }

    if (shift_period > 0) initialize_shifting(motif_models);

    ofstream sites_output_file("independent_walks.txt");
//    FILE* sites_output_file = fopen("independent_walks.txt", "w+");
    //FILE* sites_output_file = stderr;

    for (independent_walk = 0; independent_walk < total_independent_walks; independent_walk++) {
	    fprintf(stdout, "doing walk: %d\n", independent_walk);

        if (seed < 0) {
            seed = time(NULL);
            fprintf(stderr, "seeding: %d\n", seed);
            dsfmt_gv_init_gen_rand(seed);
        } else {
            dsfmt_gv_init_gen_rand(seed + iteration + independent_walk);
            fprintf(stderr, "seeding: %d\n", seed + iteration + independent_walk);
        }

        calculate_site_counts(sequences, motif_models, blocks, max_sites);

        if (starting_from_checkpoint != 1) {
            starting_from_checkpoint = 0;
            if (!sites_from_arguments) {
                cerr << "sampling initial sites uniform random" << endl;
                /**
                 * select initial random samples
                 */
                for (unsigned int i = 0; i < sequences.size(); i++) {
                    sequences.at(i).sample_uniform_random(motif_models);
                }
            } else {
                cerr << "sites were from arguments" << endl;
            }
        }

        long total_number_nucleotides = 0;
        for (unsigned int i = 0; i < sequences.size(); i++) {
            sequences.at(i).zero_accumulated_samples();
            sequences.at(i).calculate_background_site_probabilities(motif_models);

            total_number_nucleotides += sequences.at(i).nucleotides.size();
        }
        cerr << sequences.size() << " sequences with " << total_number_nucleotides << " total base pairs." << endl;

        /**
         * calculate the counts for all the samples -- we don't need to do this at 
         * every step, we can subtract the left out sequence from the counts then 
         * add the newly sampled counts
         */
        for (unsigned int i = 0; i < sequences.size(); i++) {
            increment_counts(motif_models, sequences.at(i));
        }
//        fprintf(stderr, "incremented counts for [%d] sequences.\n", number_sequences);

        printf("burn in period: %d, sample period: %d\n", burn_in_period, sample_period);
        /**
         * do the burn in walk
         */

        print_current_sites_frequency = print_current_sites_frequency_initial;

        if (total_independent_walks > 1) sites_output_file << "<independent_walk>" << endl;
        if (print_current_sites_frequency > 0) {
            sites_output_file << "<current_sites>" << endl;
            sites_output_file << "<iteration>0</iteration>" << endl;
            write_sites_to_file(sites_output_file, ".\n", sequences);
            sites_output_file << "</current_sites>" << endl;
            sites_output_file.flush();
        }

        for (int i = iteration; i < burn_in_period + sample_period; i++) {

            if (i > 0 && shift_period > 0 && (i % shift_period) == 0) attempt_shifting(max_shift_distance, sequences, motif_models);

            for (unsigned int j = 0; j < sequences.size(); j++) {
                Sequence *sequence = &(sequences.at(j));
                vector<Sample> *current_samples = &(sequence->sampled_sites);

                /**
                 * leave one sequence out and re-calculate the motif models
                 * we can do this by decrementing the counts of the left-out sequence, and recalculating the models
                 */
                decrement_counts(motif_models, *sequence);
               
                /**
                 * update the models
                 */
                update_motif_models(motif_models);

                /**
                 * resample within the sequence left out
                 */
                sequence->calculate_site_probabilities(motif_models, max_sites);

                sequence->resample_from_models(motif_models);
                
                /**
                 * update the counts with the new sample
                 */
                increment_counts(motif_models, *sequence);

                /**
                 * update the models
                 */
                update_motif_models(motif_models);

                if (i >= burn_in_period) {
                    //add the samples taken this iteration to the saved samples

                    for (unsigned int k = 0; k < current_samples->size(); k++) {
                        if (current_samples->at(k).end_position < 0 || current_samples->at(k).motif_model < 0) continue;
                        sequence->accumulated_samples.at(current_samples->at(k).motif_model).at(current_samples->at(k).end_position)++;
                    }
                }
            
                progress = ((double)i)/((double)(burn_in_period + sample_period));
                progress += (((double)j)/((double)number_sequences) * (1.0/((double)(burn_in_period + sample_period))));

#ifdef _BOINC_
                boinc_fraction_done(progress);
#else
//                fprintf(stdout, "\r%lf", progress);
#endif
            }

            if (i % 5000 == 0 && i != 0) {
                write_sites(string(SITES_CHECKPOINT_FILE), sequences, seed, i + 1, independent_walk);
                
                if (i >= burn_in_period) write_accumulated_samples(string(SAMPLES_CHECKPOINT_FILE), sequences);

                dsfmt_gv_init_gen_rand(seed + i + 1);
#ifdef _BOINC_
                boinc_checkpoint_completed();
#endif
            }

            if (print_current_sites_frequency > 0 && (((i+1) % print_current_sites_frequency) == 0)) {
                sites_output_file << "<current_sites>" << endl;
                sites_output_file << "<iteration>" << i + 1 << "</iteration>" << endl;
                write_sites_to_file(sites_output_file, ".\n", sequences);
                sites_output_file << "</current_sites>" << endl;
                sites_output_file.flush();

                if (print_current_sites_logarithmic) print_current_sites_frequency *= 2;
            }
        }

        printf("\n");
    //        printf("samples for sequence [%s]\n", sequences[i]->nucleotides);

        if (print_best_sites > 0) {
            for (unsigned int j = 0; j < motif_models.size(); j++) {
                for (unsigned int i = 0; i < sequences.size(); i++) {

                    for (unsigned int k = 0; k < sequences.at(i).nucleotides.size(); k++) {
                        count_percentage = ((double)sequences.at(i).accumulated_samples.at(j).at(k)) / ((double)sample_period);

                        if (count_percentage > print_best_sites) {
                            printf("%5d,%2d%8d ", i, j, (k + 1) - motif_models.at(j).motif_width);
                            print_sample_and_nearest(sequences.at(i), motif_models.at(j), k, 5);
                            printf(" %7d %2.2lf %s\n", k, count_percentage, sequences.at(i).name.c_str());
                        }
                    }
                }
                printf("\n");
            }
        }

        if (print_current_sites) {
            fprintf(stderr, "<current_sites>\n");
            write_sites_to_file(cerr, ".\n", sequences);
            fprintf(stderr, "</current_sites>\n");
        }

        if (print_accumulated_samples) {
            fprintf(stderr, "<current_samples>\n");
            write_accumulated_samples_to_file(cerr, sequences);
            fprintf(stderr, "</current_samples>\n");
        }

        if (sample_period > 0) write_accumulated_samples(string(SAMPLES_OUTPUT_FILE), sequences);

        if (print_motif_models) {
            for (unsigned int i = 0; i < motif_models.size(); i++) {
                motif_models.at(i).print(cerr);
            }
        }

        if (total_independent_walks > 1) sites_output_file << "</independent_walk>" << endl << endl;
    }
    sites_output_file.close();

    #ifdef _BOINC_
        boinc_finish(0);
    #endif

    return 0;
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode){
    LPSTR command_line;
    char* argv[100];
    int argc;

    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}
#endif
