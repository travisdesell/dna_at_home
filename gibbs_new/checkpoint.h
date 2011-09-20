#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include <stdio.h>

void write_sites(const char *filename, int seed, int iteration);

void read_sites_from_string(char* sites_string);

int read_sites(const char *filename);
int read_sites_from_checkpoint(const char *filename, int *seed, int *iteration);

void write_accumulated_samples_to_file(FILE* file);
void write_accumulated_samples(const char *filename);
int read_accumulated_samples(const char *filename);

#endif
