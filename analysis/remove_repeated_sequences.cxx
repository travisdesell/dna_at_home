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

    ifstream infile(input_file);

    //cout << "sequence name, character, number repeats, offset" << endl;
    string sequence_name, sequence, empty_line;
    while (infile.good()) {
       getline(infile, sequence_name); 
       //cout << "sequence_name is: " << sequence_name << endl;
       cout << sequence_name << endl;

       getline(infile, sequence); 
       //cout << "sequence is: " << sequence << endl;

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
    }
}
