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

#define MODEL_TYPE_NORMAL 0
#define MODEL_TYPE_REVERSE 1
#define MODEL_TYPE_FORWARD 2
#define MODEL_TYPE_PALINDROMIC 3

typedef struct {                                                                        

    int type;
    int motif_width;

    /**
     * counts is sized [motif_width][ALPHABET_LENGTH], it contains the sum of each letter found by the motif model at its sampled positions
     */
    double **counts;
                                                                                                         
    /**
     * nucleotide_probabilities is sized [motif_width][ALPHABET_LENGTH], it contains the probability of finding a letter according to counts
     */
    double **nucleotide_probabilities;
    double **reverse_nucleotide_probabilities;
} MotifModel;                                                                                            
                                                                                                         

typedef struct {
    int end_position;
    int motif_model;
} Sample;

typedef struct {
    char *name;

    int length;
    char *nucleotides;
    int *distance_to_previous_invalid;      // used for calculating if a position is a possible end position for a motif model

    long double **possible_site_counts;       // counts the number of possible ways to place [1 .. max_sites] motif model sites within the sequence

    long double *possible_sites;                //the values of the last column of possible site counts
    long double total_possible_sites;           //the sum of possible_sites

    long double **background_site_probability;
    long double **site_probability_ratio;

    long double **site_probability;
    long double ***motif_probability_contribution;

    Sample **sampled_sites;
    Sample **shifted_sites;

} Sequence;

typedef struct {
    char *name;
    double evolution_time;
    long double *nucleotide_probability;
    long double **matrix;

} PhylogenyTreeNode;


#endif
