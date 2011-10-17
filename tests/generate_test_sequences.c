#include <stdlib.h>
#include <time.h>

#include "../../mersenne_twister/dSFMT.h"


int main(int number_arguments, char **arguments) {
    int character, i, j, k;
    char **sequences, **motifs;
    int **start_position;
    int number_sequences, sequence_length;
    int number_motifs, motif_length;
    int seed;

    printf("arguments[1]: %s\n", arguments[1]);

    number_sequences = -1;
    sequence_length = -1;
    number_motifs = -1;
    motif_length = -1;

    seed = time(NULL);

    for (i = 1; i < number_arguments; i++) {
        if ( !strcmp(arguments[i], "--number_sequences") ) {
            number_sequences = atoi(arguments[++i]);

        } else if ( !strcmp(arguments[i], "--sequence_length") ) {
            sequence_length = atoi(arguments[++i]);

        } else if ( !strcmp(arguments[i], "--number_motifs") ) {
            number_motifs = atoi(arguments[++i]);

        } else if ( !strcmp(arguments[i], "--motif_length")) {
            motif_length = atoi(arguments[++i]);

        } else if ( !strcmp(arguments[i], "--seed")) {
            seed = atoi(arguments[++i]);

        }
    }


    sequences = (char**)malloc(number_sequences * sizeof(char*));
    for (i = 0; i < number_sequences; i++) sequences[i] = (char*)malloc((sequence_length + 1) * sizeof(char));

    motifs = (char**)malloc(number_motifs * sizeof(char*));
    for (i = 0; i < number_motifs; i++) motifs[i] = (char*)malloc((motif_length + 1) * sizeof(char));

    start_position = (int**)malloc(number_sequences * sizeof(int*));
    for (i = 0; i < number_sequences; i++) {
        start_position[i] = (int*)malloc(number_motifs * sizeof(int));
        for (j = 0; j < number_motifs; j++) {
            start_position[i][j] = -1;
        }
    }

    dsfmt_gv_init_gen_rand(seed);

    for (i = 0; i < number_sequences; i ++) {
        for (j = 0; j < sequence_length; j++) {
            character = (int)(dsfmt_gv_genrand_open_close() * 4.0);

            if (character == 0) sequences[i][j] = 'A';
            if (character == 1) sequences[i][j] = 'C';
            if (character == 2) sequences[i][j] = 'G';
            if (character == 3) sequences[i][j] = 'T';
        }

        sequences[i][sequence_length] = '\0';
//        printf(">sequence[%d]\n", i);
//        printf("%s\n\n", sequences[i]);
    }

    for (i = 0; i < number_motifs; i++) {
        for (j = 0; j < motif_length; j++) {
            character = (int)(dsfmt_gv_genrand_open_close() * 4.0);

            if (character == 0) motifs[i][j] = 'A';
            if (character == 1) motifs[i][j] = 'C';
            if (character == 2) motifs[i][j] = 'G';
            if (character == 3) motifs[i][j] = 'T';
        }
        motifs[i][motif_length] = '\0';
        printf("motif [%d]: [%s]\n", i, motifs[i]);
    }
//    printf("overrepresented sequence: [%s]\n", motifs);

    for (i = 0; i < number_sequences; i++) {
        for (j = 0; j < number_motifs; j++) {
            start_position[i][j] = (int)(dsfmt_gv_genrand_open_close() * ((sequence_length - motif_length) - (j* motif_length)));
//        printf("positions in sequence [%d]: %d to %d\n", i, start_position[i], start_position[i] + motif_length);

            if (j > 0 && start_position[i][j-1]) {
                if ((start_position[i][j] >= start_position[i][j-1]) && (start_position[i][j] < start_position[i][j-1] + motif_length)) {
                    printf("[%d] start position %d collided with start position %d\n", i, start_position[i][j], start_position[i][j-1]);
                    start_position[i][j] += motif_length;
                }
            }

            for (k = 0; k < motif_length; k++) {
                sequences[i][start_position[i][j] + k] = motifs[j][k];
            }
        }
    }

//    printf("modified sequences:\n");
    for (i = 0; i < number_sequences; i++) {
        printf(">sequence[%d]", i);
        for (j = 0; j < number_motifs; j++) {
            printf(", start pos [%d] subseq [%d][%s] end pos [%d]", start_position[i][j], j, motifs[j], start_position[i][j] + motif_length - 1);
        }
        printf("\n");

        printf("%s\n\n", sequences[i]);
    }
}
