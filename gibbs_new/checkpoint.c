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

#include "structs.h"
#include "checkpoint.h"
#include "sequences.h"


extern int max_sites;

extern int number_sequences;
extern Sequence **sequences;

extern int number_motifs;

extern int ***accumulated_samples;

void write_sites(const char *filename, int seed, int iteration) {
    FILE *checkpoint_file;
    int retval;

#ifdef _BOINC_
    char output_path[512];
    boinc_resolve_filename(filename, output_path, sizeof(output_path));

    checkpoint_file = boinc_fopen(output_path, "w+");
#else
    checkpoint_file = fopen(filename, "w+");
#endif
    if (!checkpoint_file) {
        fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file)\n");
        return;
    }

    fprintf(checkpoint_file, "seed: %d\n", seed);
    fprintf(checkpoint_file, "iteration: %d\n", iteration);
    write_sites_to_file(checkpoint_file, ".\n");

    if ((retval = fclose(checkpoint_file))) {
        fprintf(stderr, "APP: error writing checkpoint (closing checkpoint file) %d\n", retval);
        return;
    }
}

void read_sites_from_string(char *sites_string) {
    int i, j, k;
    int current_sequence, current_sample;
    int sites_string_length;
    char current_sequence_string[50];
    char *token;
    Sample *sample;

    fprintf(stderr, "reading sites from string: [%s]\n", sites_string);

    sites_string_length = strlen(sites_string);
    printf("sites string length: %d\n", sites_string_length);

    current_sequence = 0;
    i = 0;
    while (i < sites_string_length) {
        if (sites_string[i] == ';') {
            current_sequence++;
            i++;

            if (i < sites_string_length) {
                for (k = 0; k < max_sites; k++) {
                    sequences[current_sequence]->sampled_sites[k]->end_position = -1;
                    sequences[current_sequence]->sampled_sites[k]->motif_model = -1;
                }
            }

            fprintf(stderr, "incrementing sequence to : %d\n", current_sequence);
            continue;
        }

        j = 0;
        while (sites_string[i + j] != ';') {
            current_sequence_string[j] = sites_string[i+j];
            j++;
        }
        current_sequence_string[j] = '\0';

        printf("current sequence string: [%s]\n", current_sequence_string);

        token = strtok(current_sequence_string, ":");

        current_sample = 0;
        while (token != NULL) {
            sample = sequences[current_sequence]->sampled_sites[current_sample];

            printf("sscanfing on token: [%s]\n", token);

            sscanf(token, "%d,%d", &(sample->motif_model), &(sample->end_position));

            token = strtok(NULL, ":\0");
            current_sample++;
        }

        for (k = 0; k < 50; k++) current_sequence_string[k] = '\0';
        i += j;
    }
}

void read_sites_from_file(FILE* file) {
    int i;
    int current_sequence, current_sample;
    int input_motif_model, input_end_position;
    char line[200];
    char *token;

    current_sequence = 0;
//    fprintf(stderr, "current_sequence: %d, number_sequences: %d\n", current_sequence, number_sequences);

    while (NULL != fgets(line, 200, file)) {
//        fprintf(stdout, "READ SITES FILE LINE: [%s]\n", line);

        if (current_sequence >= number_sequences) {
            fprintf(stderr, "ERROR reading sites from file. Read more sequences than number_sequences, In [%s], line [%d]\n", __FILE__, __LINE__);
            fprintf(stderr, "current_sequence: %d, number_sequences: %d\n", current_sequence, number_sequences);
            exit(0);
        }

        for (i = 0; i < max_sites; i++) {
            sequences[current_sequence]->sampled_sites[i]->end_position = -1;
            sequences[current_sequence]->sampled_sites[i]->motif_model = -1;
        }

        current_sample = 0;
        token = strtok(line, ".,:\r\n");
        while (NULL != token) { 
            if (current_sample >= max_sites) {
                fprintf(stderr, "ERROR reading sites from file.  Read more sites than max_sites. In [%s], line [%d]\n", __FILE__, __LINE__);
                fprintf(stderr, "current_sample: %d, max_sites: %d\n", current_sample, max_sites);
                fprintf(stderr, "line from file is: [%s]\n", line);
                fprintf(stderr, "token is: [%s]\n", token);

                exit(0);
            }

//            printf("setting sequence[%d]->sampled_sites[%d]\n", current_sequence, current_sample);

//            fprintf(stdout, "token: [%s], int value: [%d]\n", token, (int)token[0]);
            input_motif_model = atoi(token);
            if (input_motif_model < 0 || input_motif_model >= number_motifs) {
                fprintf(stderr, "ERROR reading sites from file.  input motif model > number motifs. In [%s], line [%d]\n", __FILE__, __LINE__);
                fprintf(stderr, "input_motif_model: %d, number_motifs: %d\n", input_motif_model, number_motifs);
                fprintf(stderr, "line from file is: [%s]\n", line);
                fprintf(stderr, "token is: [%s]\n", token);

                exit(0);
             }

            sequences[current_sequence]->sampled_sites[current_sample]->motif_model = input_motif_model;
            token = strtok(NULL, ",:\r\n");

//            fprintf(stdout, "token: [%s], int value: [%d]\n", token, (int)token[0]);
            input_end_position = atoi(token);
            if (input_end_position < 0 || input_end_position >= sequences[current_sequence]->length) {
                fprintf(stderr, "ERROR reading sites from file.  input end position > current_sequence->length. In [%s], line [%d]\n", __FILE__, __LINE__);
                fprintf(stderr, "input_end_position: %d, current_sequence->length: %d\n", input_end_position, sequences[current_sequence]->length);
                fprintf(stderr, "line from file is: [%s]\n", line);
                fprintf(stderr, "token is: [%s]\n", token);

                exit(0);
             }


            sequences[current_sequence]->sampled_sites[current_sample]->end_position = input_end_position;
            token = strtok(NULL, ",:\r\n");
            current_sample++;
        }

        current_sequence++;
    }

    if (current_sequence != number_sequences) {
            fprintf(stderr, "ERROR reading sites from file. Did not read sites for each sequence. In [%s], line [%d]\n", __FILE__, __LINE__);
            fprintf(stderr, "current_sequence: %d, number_sequences: %d\n", current_sequence, number_sequences);
            exit(0);
    }
}


int read_sites(const char *filename) {
    FILE *file;

#ifdef _BOINC_    
    char input_path[512];
    int retval;

    retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
    if (retval) {
        return 0;
    }

    file = boinc_fopen(input_path, "r");
#else
    file = fopen(filename, "r");
#endif

    if (file == NULL) {
        printf("sites file was null\n");
        //No checkpoint was found
        return 0;
    }

    printf("reading from sites from arguments\n");

    read_sites_from_file(file);

    fclose(file);

    return 1;
}

int read_sites_from_checkpoint(const char *filename, int *seed, int *iteration) {
    FILE *checkpoint_file;

#ifdef _BOINC_    
    char input_path[512];
    int retval;

    retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
    if (retval) {
        return 0;
    }

    checkpoint_file = boinc_fopen(input_path, "r");
#else
    checkpoint_file = fopen(filename, "r");
#endif

    if (checkpoint_file == NULL) {
//        printf("sites checkpoint file was null\n");
        //No checkpoint was found
        return 0;
    }

    printf("reading from sites checkpoint\n");

    fscanf(checkpoint_file, "seed: %d\n", seed);
    fscanf(checkpoint_file, "iteration: %d\n", iteration);

    read_sites_from_file(checkpoint_file);

    fclose(checkpoint_file);

    return 1;
}


void write_accumulated_samples_to_file(FILE *file) {
    int i, j, k;

    for (i = 0; i < number_sequences; i++) {
        for (j = 0; j < number_motifs; j++) {
            fprintf(file, "%d", accumulated_samples[i][j][0]);
            for (k = 1; k < sequences[i]->length; k++) {
                fprintf(file, ",%d", accumulated_samples[i][j][k]);
            }
            fprintf(file, "\n");
        }
    }
}

void write_accumulated_samples(const char *filename) {
    FILE *checkpoint_file;

#ifdef _BOINC_
    char output_path[512];

    boinc_resolve_filename(filename, output_path, sizeof(output_path));

    checkpoint_file = boinc_fopen(output_path, "w+");
    if (!checkpoint_file) {
        fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file)\n");
        return;
    }
#else
    checkpoint_file = fopen(filename, "w+");
#endif

    write_accumulated_samples_to_file(checkpoint_file);
    fclose(checkpoint_file);
}

int read_accumulated_samples(const char *filename) {
    int i, j, k;
    FILE *checkpoint_file;
    char line[10000];
    char *token;

#ifdef _BOINC_    
    char input_path[512];
    int retval;

    retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
    if (retval) {
        return 0;
    }

    checkpoint_file = boinc_fopen(input_path, "r");
#else
    checkpoint_file = fopen(filename, "r");
#endif

    if (checkpoint_file == NULL) {
        fprintf(stderr, "samples checkpoint file was null\n");
        //No checkpoint was found
        return 0;
    }
    fprintf(stderr, "reading from samples checkpoint\n");

    for (i = 0; i < number_sequences; i++) {
        for (j = 0; j < number_motifs; j++) {
            fgets(line, 10000, checkpoint_file);

            if (line == NULL) {
                fprintf(stderr, "ERROR in reading samples checkpoint, should read more sequence information but line == NULL\n");
                fprintf(stderr, "ERROR on line [%d], file [%s]\n", __LINE__, __FILE__);
                exit(0);
            }

//            fprintf(stderr, "READ LINE: [%s] length [%d]\n", line, (int)strlen(line));
//            fprintf(stderr, "sequences[%d]->length: [%d]\n", i, sequences[i]->length);

            token = strtok(line, ",\r\n");

            for (k = 0; k < sequences[i]->length; k++) {
                if (token == NULL) {
                    fprintf(stderr, "ERROR: reading samples, token == NULL before all samples should have been read.\n");
                    fprintf(stderr, "loop i: [%d], j: [%d], k: [%d]\n", i, j, k); 
                    fprintf(stderr, "error on line [%d], file [%s]\n", __LINE__, __FILE__);
                    exit(0);
                }
//                printf("token: [%s], k: [%d]\n", token, k);

                accumulated_samples[i][j][k] = atoi(token);

                token = strtok(NULL, ",\r\n");
            }
            if (token != NULL) {
                fprintf(stderr, "ERROR: reading samples, token != NULL after all samples should have been read.\n");
                fprintf(stderr, "token is: [%s]\n", token);
                fprintf(stderr, "next token is: [%s]\n", strtok(NULL, ",\r\n"));
                fprintf(stderr, "i: [%d], j: [%d]\n", i, j);
                fprintf(stderr, "error on line [%d], file [%s]\n", __LINE__, __FILE__);
                exit(0);
            }
        }
    }

    fclose(checkpoint_file);
    return 1;
}
