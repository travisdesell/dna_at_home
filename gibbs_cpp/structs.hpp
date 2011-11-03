#include <iostream>
#include <string>
#include <vector>

using namespace std;

#ifndef STRUCTS_H
#define STRUCTS_H

#define ALPHABET_LENGTH 4
#define ALPHABET "ACGT"

#define A_POSITION 0
#define C_POSITION 1
#define G_POSITION 2
#define T_POSITION 3

#define A_POSITION_REVERSE 3
#define C_POSITION_REVERSE 2
#define G_POSITION_REVERSE 1
#define T_POSITION_REVERSE 0

/**
 *  This is a prime example of why C sucks.  Can't have the structs spread out across multiple header files.
 */
class PhylogenyTreeNode {
    char *name;
    double evolution_time;
    vector<long double> nucleotide_probability;
    vector< vector<long double> > matrix;
};


#endif
