#include <cstdlib>
#include <stdint.h>

#include <fstream>
using std::ifstream;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <sstream>
using std::ostringstream;

#include <string>
using std::string;


int main(int n_arguments, char **arguments) {
    if (n_arguments != 3) {
        cerr << "ERROR: incorrect arguments. Usage:" << endl;
        cerr << "\t" << arguments[0] << " <input file (FASTA)> <min repeat length>" << endl;
        exit(1);
    }

    string input_file(arguments[1]);
    int min_repeat_length = atoi(arguments[2]);

    ifstream infile(input_file.c_str());

    //cout << "sequence name, character, number repeats, offset" << endl;
    string sequence_name, sequence, empty_line;
    getline(infile, sequence_name); 
    while (infile.good()) {
        if (sequence_name.size() > 0 && sequence_name[0] != '>') {
            cout << "invalid sequence name: '" << sequence_name << "'" << endl;
            break;
        }
        cout << sequence_name << endl;

        if (infile.bad()) break;
        getline(infile, sequence); 
        if (sequence.size() > 0
                && sequence[0] != 'A'
                && sequence[0] != 'C'
                && sequence[0] != 'G'
                && sequence[0] != 'T'
                && sequence[0] != 'a'
                && sequence[0] != 'c'
                && sequence[0] != 'g'
                && sequence[0] != 't'
                && sequence[0] != 'X'
                && sequence[0] != 'x') {
            cout << "invalid sequence: " << sequence << endl;
            break;
        }

        getline(infile, empty_line);

        char repeated_character = toupper(sequence[0]);
        int repeat_count = 1;
        int repeat_start = 0;

        ostringstream oss;

        for (uint32_t i = 1; i < sequence.size(); i++) {
            char current_character = toupper(sequence[i]);
            if (current_character == repeated_character) {
                repeat_count++;
            } else {
                if (repeat_count >= min_repeat_length) {
                    //cout << "repeated sequence of " << repeat_count << " " << repeated_character << " started at " << repeat_start << ", offset from gene start: " << (repeat_start - offset) << endl;
                    //cout << sequence_name.substr(2, sequence_name.size() - 2) << ", " << current_character << ", " << repeat_count << ", "<< (repeat_start - offset) << endl;
                    for (uint32_t j = 0; j < repeat_count; j++) {
                        oss << "X";
                    }
                } else {
                    for (uint32_t j = 0; j < repeat_count; j++) {
                        oss << repeated_character;
                    }
                }

                repeated_character = current_character;
                repeat_count = 1;
                repeat_start = i;
            }
        }
        for (uint32_t j = 0; j < repeat_count; j++) {
            oss << repeated_character;
        }


        cout << oss.str() << endl;
        if (oss.str().size() != sequence.size()) {
            cout << "sequence.size(): " << sequence.size() << " != oss.str().size(): " << oss.str().size() << endl;
            cout << "sequence:" << endl;
            cout << sequence << endl;
            exit(1);
        }
        cout << endl;

        getline(infile, sequence_name); 
    }
}
