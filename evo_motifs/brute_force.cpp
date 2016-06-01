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

#include <unordered_map>
using std::unordered_map;

#include <vector>
using std::vector;

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

string reverse_complement(string motif) {
    string rc; 

    for (int32_t i = motif.size() - 1; i >= 0; i--) {
        if (motif[i] == 'A') rc += 'T';
        else if (motif[i] == 'C') rc += 'G';
        else if (motif[i] == 'G') rc += 'C';
        else if (motif[i] == 'T') rc += 'A';
        else if (motif[i] == 'X') rc += 'X';
    }   

    //cout << "motif '" << motif << "' -- reverse complement '" << rc << "'" << endl;

    return rc; 
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

class Motif {
    private:
        int forward_count;
        int reverse_count;

        string forward;
        string reverse;

    public:

        int get_forward_count() const { return forward_count; }
        int get_reverse_count() const { return reverse_count; }
        int get_total_count() const { return forward_count + reverse_count; }

        int distinct_nucleotide_count() const {
            int A_count = 0;
            int C_count = 0;
            int G_count = 0;
            int T_count = 0;

            for (uint32_t i = 0; i < forward.size(); i++) {
                if (forward[i] == 'A') A_count = 1;
                else if (forward[i] == 'C') C_count = 1;
                else if (forward[i] == 'G') G_count = 1;
                else if (forward[i] == 'T') T_count = 1;
            }

            return A_count + C_count + G_count + T_count;
        }

        Motif() {
            forward_count = 0;
            reverse_count = 0;
            forward = "";
            reverse = "";
        }

        Motif(string nucleotides) {
            forward_count = 1;
            reverse_count = 0;

            forward = nucleotides;
            reverse = reverse_complement(nucleotides);
        }

        void increment(string nucleotides) {
            if (nucleotides.compare(forward) == 0) {
                forward_count++;
            } else if (nucleotides.compare(reverse) == 0) {
                reverse_count++;
            } else {
                cerr << "ERROR! incremeting a motif:" << endl;
                cerr << "\tforward: '" << forward << "' : " << forward_count << endl;
                cerr << "\treverse: '" << reverse << "' : " << reverse_count << endl;
                cerr << "with string '" << nucleotides << "' that does not match forward or reverse!" << endl;
                exit(1);
            }
        }

        friend ostream &operator<<(std::ostream &os, Motif const &motif);
};

ostream &operator<<(std::ostream &os, Motif const &motif) {
//    os << "total: " << motif.get_total_count() << ", forward: '" << motif.forward << "' : " << motif.forward_count << ", reverse: '" << motif.reverse << "' : " << motif.reverse_count;
    os << motif.get_total_count() << ", " << motif.forward << ", " << motif.forward_count << ", " << motif.reverse << ", " << motif.reverse_count;
    return os;
}

Sequences sequences;

struct sort_motifs_desc {
    bool operator()(const Motif &motif1, const Motif &motif2) {
        return motif1.get_total_count() > motif2.get_total_count();
    }   
};


int main(int argc, char **argv) {
    vector<string> arguments(argv, argv + argc);

    string sequences_filename = arguments[1];
    string position_weight_matrix_filename = arguments[2];

    sequences.read_from_file(sequences_filename);
    sequences.set_verbose(false);

    //cout << sequences << endl;

    uint32_t motif_size = 6;
    unordered_map<string, Motif> motif_counts;

    for (uint32_t i = 0; i < sequences.number_sequences(); i++) {
        string nucleotides = sequences.get(i).get_nucleotides();

        //cout << "NUCLEOTIDES: " << nucleotides << endl;

        string current = nucleotides.substr(0, motif_size);

        if (motif_counts.count(current) > 0) {
            motif_counts[current].increment(current);
        } else if (motif_counts.count( reverse_complement(current) ) > 0) {
            motif_counts[reverse_complement(current)].increment(current);
        } else {
            motif_counts.insert({current, Motif(current)});
        }
        //cout << "current: " << motif_counts[current] << endl;

        for (uint32_t j = motif_size; j < nucleotides.size(); j++) {
            current.erase(0, 1);
            current += nucleotides[j];

            if (motif_counts.count(current) > 0) {
                motif_counts[current].increment(current);
            } else if (motif_counts.count( reverse_complement(current) ) > 0) {
                motif_counts[reverse_complement(current)].increment(current);
            } else {
                motif_counts.insert({current, Motif(current)});
            }
            //cout << "current: " << motif_counts[current] << endl;
        }

        //cout << endl << endl;
    }

    //cout << "MOTIF LIST" << endl;
    vector<Motif> motifs;
    motifs.reserve(motif_counts.size());

    for(auto kv : motif_counts) {
        motifs.push_back(kv.second);  
        //cout << kv.second << endl;
    } 

    cout << "#SORTED MOTIF LIST" << endl;
    cout << "#total count, foward motif, forward count, reverse motif, reverse count" << endl;
    sort(motifs.begin(), motifs.end(), sort_motifs_desc());

    for (uint32_t i = 0; i < motifs.size(); i++) {
        if (motifs[i].get_total_count() > 1 && motifs[i].distinct_nucleotide_count() > 3) {
            cout << motifs[i] << endl;
        }
    }


}
