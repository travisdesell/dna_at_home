#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef _BOINC_
    #include "diagnostics.h"
    #include "util.h"
    #include "filesys.h"
    #include "boinc_api.h"
    #include "mfile.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "structs.hpp"
#include "checkpoint.hpp"
#include "sequences.hpp"

void write_sites(string filename, const vector<Sequence> &sequences, int seed, int iteration, int independent_walk) {
#ifdef _BOINC_
    string output_path;
    int retval = boinc_resolve_filename_s(filename, output_path);
    
    if (retval) {
        fprintf(stderr, "APP: error writing checkpoint (resolving checkpoint file name)\n");
        return;
    }

    ofstream checkpoint_file(output_path.c_str());
#else
    ofstream checkpoint_file(filename.c_str());
#endif
    if (!checkpoint_file.is_open()) {
        fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file)\n");
        return;
    }

    checkpoint_file << "seed: " <<  seed << endl;
    checkpoint_file << "iteration: " << iteration << endl;
    checkpoint_file << "independent_walk: " << independent_walk << endl;
    write_sites_to_file(checkpoint_file, ".\n", sequences);

    checkpoint_file.close();
}

void write_sites_to_file(ostream &out, const char *delimiter, const vector<Sequence> &sequences) {
    bool already_printed = false;

    for (unsigned int i = 0; i < sequences.size(); i++) {
        for (unsigned int k = 0; k < sequences.at(i).sampled_sites.size(); k++) {
            const Sample *sample = &(sequences.at(i).sampled_sites.at(k));
            if (sample->end_position >= 0) {
                if (already_printed) out << ":";
                out << sample->motif_model << "," << sample->end_position;
                already_printed = true;
            }   
        }   
        out << delimiter;
        already_printed = false;
    }   
}


void read_sites(string sites_filename, vector<Sequence> &sequences) {
#ifdef _BOINC_
    string input_path;
    int retval;

    retval = boinc_resolve_filename_s(filename, input_path);
    if (retval) {
        return 0;
    }

    ifstream sites_file(input_path);
#else 
    ifstream sites_file(sites_filename.c_str());
#endif

    read_sites_from_file(sites_file, sequences);
}

int read_sites_from_checkpoint(string sites_filename, vector<Sequence> &sequences, int &seed, int &iteration, int &independent_walk) {
#ifdef _BOINC_
    char input_path[512];
    int retval;

    retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
    if (retval) {
        return 0;
    }

    ifstream sites_file(input_path); /** TODO: NEED TO DOUBLE CHECK IF THIS IS POSSIBLE. **/
#else 
    ifstream sites_file(sites_filename.c_str());
#endif
    if (!sites_file.is_open()) return 0;

    string s;
    sites_file >> s >> seed;
    if (s.compare("seed:") != 0) {
        cerr << "ERROR: malformed sites checkpoint! could not read 'seed'" << endl;
        exit(0);
    }
//    cerr << "seed: " << seed << endl;

    sites_file >> s >> iteration;
    if (s.compare("iteration:") != 0) {
        cerr << "ERROR: malformed sites checkpoint! could not read 'iteration'" << endl;
        exit(0);
    }
//    cerr << "iteration: " << iteration << endl;

    sites_file >> s >> independent_walk;
    if (s.compare("independent_walk:") != 0) {
        cerr << "ERROR: malformed sites checkpoint! could not read 'independent_walk'" << endl;
        exit(0);
    }
//    cerr << "independent_walk: " << independent_walk << endl;

    getline(sites_file, s); //need this because the >> does not grab the newline

    read_sites_from_file(sites_file, sequences);
    return 1;
}

void read_sites_from_file(ifstream &sites_file, vector<Sequence> &sequences) {
    if (sites_file.is_open()) {
        string current_line;
        int current_sequence = 0;

        getline(sites_file, current_line);
        while (sites_file.good() && !sites_file.eof()) {
//            cout << "current sequence: " << current_sequence << ", current line: " << current_line << endl;
            if (current_sequence >= (int)sequences.size()) {
                cerr << "ERROR: sites_file contains sites for more sequences than in sequences file." << endl;
                exit(0);
            }

            sequences.at(current_sequence).sampled_sites.clear();

            int motif_number, end_position;
            stringstream ss(current_line);
            char c;
            while ( ss >> motif_number >> c  >> end_position ) {
                ss >> c;

                sequences.at(current_sequence).sampled_sites.push_back(Sample(motif_number, end_position));

                if (c == '.') break;
            }
            current_sequence++;
            getline(sites_file, current_line);
        }
        sites_file.close();
    }
}

void write_accumulated_samples_to_file(ostream &out, const vector<Sequence> &sequences) {
    for (unsigned int i = 0; i < sequences.size(); i++) {
        const Sequence *sequence = &(sequences.at(i));
        for (unsigned int j = 0; j < sequence->accumulated_samples.size(); j++) {
            out << sequence->accumulated_samples.at(j).at(0);
            for (unsigned int k = 1; k < sequence->accumulated_samples.at(j).size(); k++) {
                out << "," << sequence->accumulated_samples.at(j).at(k);
            }
            out << endl;
        }
    }
}

void write_accumulated_samples(string filename, const vector<Sequence> &sequences) {
#ifdef _BOINC_
    string output_path;

    retval = boinc_resolve_filename_s(filename, output_path);
    if (!checkpoint_file) {
        fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file)\n");
        return;
    }

    ofstream out(output_path.c_str());
#else
    ofstream out(filename.c_str());
#endif

    write_accumulated_samples_to_file(out, sequences);
    out.close();
}

int read_accumulated_samples(string filename, vector<Sequence> &sequences) {
#ifdef _BOINC_    
    string input_path;
    int retval;

    retval = boinc_resolve_filename_s(filename, input_path);
    if (retval) {
        return 0;
    }

    ifstream instream(filename.c_str());
#else
    ifstream instream(filename.c_str());
#endif

    if (!instream.is_open()) {
        cerr << "ERROR: could not open accumulated samples checkpoint file for reading" << endl;
        //No checkpoint was found
        return 0;
    }
    cerr << "reading from samples checkpoint" << endl;

    string line;

    for (unsigned int i = 0; i < sequences.size(); i++) {
        Sequence *sequence = &(sequences.at(i));

        for (unsigned int j = 0; j < sequence->accumulated_samples.size(); j++) {
            if ( !getline(instream, line) ) {
                cerr << "ERROR in reading samples checkpoint, should read more sequence information but could not read another line" << endl;
                cerr << "ERROR on line [" << __LINE__ << "], file [" <<__FILE__ << "]" << endl;
                exit(0);
            }

            char *line_c = (char*)malloc(line.size() + 1);
            copy(line.begin(), line.end(), line_c);
            line_c[line.size()] = '\0';

//            cout << "current line: " << line_c << endl;
            char *token = strtok(line_c, ",\r\n");

//            cout << "acc samples.at(j).size: " << sequence->accumulated_samples.at(j).size() << endl;
//            cout << "sequence nucleotides size: " << sequence->nucleotides.size() << endl;

            for (unsigned int k = 0; k < sequence->accumulated_samples.at(j).size(); k++) {
                if (token == NULL) {
                    cerr << "ERROR: reading samples, token == NULL before all samples should have been read." << endl;
                    cerr << "loop i: [" << i << "], j: [" << j << "], k: [" << k << "]" << endl;
                    cerr << "error on line [" << __LINE__ << "], file [" << __FILE__ << "]" << endl;
                    exit(0);
                }
//                printf("token: [%s], k: [%d]\n", token, k);

                sequence->accumulated_samples.at(j).at(k) = atoi(token);

                token = strtok(NULL, ",\r\n");
            }
            if (token != NULL) {
                cerr << "ERROR: reading samples, token != NULL after all samples should have been read." << endl;
                cerr << "token is: [" << token << "]" << endl;
                cerr << "next token is: [" << strtok(NULL, ",\r\n") << "]" << endl;
                cerr << "i: [" << i << "], j: [" << j << "]" << endl;
                cerr << "error on line [" << __LINE__ << "], file [" << __FILE__ << "]" << endl;
                exit(0);
            }
            delete [] line_c;
        }
    }

    instream.close();
    return 1;
}
