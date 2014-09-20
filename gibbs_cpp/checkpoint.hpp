#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include <vector>
#include <iostream>
#include <fstream>

#include "stdio.h"
#include "sequences.hpp"

using namespace std;

void write_sites(string filename, vector<Sequence*> *sequences, int seed, int iteration, int indendent_walk);
void write_final_sites(string filename, vector<Sequence*> *sequences);
void write_sites_to_file(ostream &out, const char *delimiter, vector<Sequence*> *sequences);

int read_sites_from_checkpoint(string checkpoint_filename, vector<Sequence*> *sequences, int &seed, int &iteration, int &independent_walk);
void read_sites(string sites_filename, vector<Sequence*> *sequences);
void read_sites_from_file(ifstream &sites_file, vector<Sequence*> *sequences);

void write_accumulated_samples_to_file(ostream &out, vector<Sequence*> *sequences);
void write_accumulated_samples(string filename, vector<Sequence*> *sequences);
int read_accumulated_samples(string filename, vector<Sequence*> *sequences);

#endif
