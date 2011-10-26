#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "../../mersenne_twister/dSFMT.h"


#define MODEL_TYPE_NORMAL 0
#define MODEL_TYPE_REVERSE 1
#define MODEL_TYPE_FORWARD 2
#define MODEL_TYPE_PALINDROMIC 3

int main(int number_arguments, char **arguments) {
    int character, i, j, k;
    char **sequences, **motifs;
    int **start_position;
    int number_sequences, sequence_length,palindrome_length;
    /* motif_build 1 = exact motif; 2=probabilistic random motif
    */
    int number_motifs, motif_length,motif_build,motif_model;
    int seed;

    printf("arguments[1]: %s\n", arguments[1]);

    number_sequences = -1;
    sequence_length = -1;
    number_motifs =-1;
    motif_length = -1;
    motif_model = 1;
    motif_build = 1;
    palindrome_length=8; //length of palindromic part of motif starting at 0
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

        } else if ( !strcmp(arguments[i], "--motif_model")) {
        	motif_model = atoi(arguments[++i]);

        } else if ( !strcmp(arguments[i], "--motif_build")) {
            motif_build = atoi(arguments[++i]);

 	} else if ( !strcmp(arguments[i], "--palindrome_length")) {
            palindrome_length = atoi(arguments[++i]);

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

    if((MODEL_TYPE_REVERSE == motif_model||MODEL_TYPE_NORMAL == motif_model) && 1==motif_build)
    {
    for (i = 0; i < number_motifs; i++) {
        for (j = 0; j < motif_length; j++) {
            character = (int)(dsfmt_gv_genrand_open_close() * 4.0);

            if (character == 0) motifs[i][j] = 'A';
            if (character == 1) motifs[i][j] = 'C';
            if (character == 2)	motifs[i][j] = 'G';
            if (character == 3)	motifs[i][j] = 'T';
        }
        motifs[i][motif_length] = '\0';
        printf("motif [%d]: [%s]\n", i, motifs[i]);
    }
    }
   
    //Palindromic motif
    if(MODEL_TYPE_PALINDROMIC == motif_model && 1 == motif_build)
    {
	if( 2*palindrome_length > motif_length)
             {
                printf("Palindrome lenght should be atmost half of the motif_length\n Palindrome Length=%d Motif Length=%d\n",palindrome_length, motif_length);
                return;
             }
    for (i = 0; i < number_motifs; i++) {
        for (j = 0; j < motif_length-palindrome_length; j++) {
            character = (int)(dsfmt_gv_genrand_open_close() * 4.0);

            if (character == 0)
            	{
            	motifs[i][j] = 'A';
		if(j<palindrome_length)
            	motifs[i][motif_length-j-1] = 'T';
            	}
            if (character == 1)
		{
            	motifs[i][j] = 'C';
		if(j<palindrome_length)
               	motifs[i][motif_length-j-1] = 'G';
		}
            if (character == 2)
		{
            	motifs[i][j] = 'G';
		if(j<palindrome_length)
                motifs[i][motif_length-j-1] = 'C';
            	}

            if (character == 3)
		{
            	motifs[i][j] = 'T';
		if(j<palindrome_length)
                motifs[i][motif_length-j-1] = 'A';
		}

        }
        motifs[i][motif_length] = '\0';
        printf("motif [%d]: [%s]\n", i, motifs[i]);
    }
	}

//    printf("overrepresented sequence: [%s]\n", motifs);

    for (i = 0; i < number_sequences; i++) {
        for (j = 0; j < number_motifs; j++) {
            start_position[i][j] = (int)(dsfmt_gv_genrand_open_close() * ((sequence_length - motif_length) - (j* motif_length)));
//        printf("positions in sequence [%d]: %d to %d\n", i, start_position[i], start_position[i] + motif_length);

            if (j > 0 /*&& start_position[i][j-1]*/) {
                if ((start_position[i][j] >= start_position[i][j-1]) && (start_position[i][j] < (start_position[i][j-1] + motif_length))) {
                    printf("[%d] start position %d collided with start position %d\n", i, start_position[i][j], start_position[i][j-1]);
                    start_position[i][j] += motif_length;
                }
                if((start_position[i][j] <= start_position[i][j-1]) && ((start_position[i][j]+ motif_length)> start_position[i][j-1] ))
                {
                	printf("[%d] start position %d collided with start position %d\n", i, start_position[i][j], start_position[i][j-1]);
                	if(start_position[i][j]-motif_length >= 0)
                	start_position[i][j] -= motif_length;
                	else
                	start_position[i][j] += 2*motif_length;
                }
            }
	   if(MODEL_TYPE_REVERSE == motif_model)
	   {
	     for (k = 0; k < motif_length; k++)
               {
	       if(motifs[j][motif_length-k-1] == 'A')
		 sequences[i][start_position[i][j] + k] = 'T';
               if(motifs[j][motif_length-k-1]=='T')
                 sequences[i][start_position[i][j] + k] = 'A';
               if(motifs[j][motif_length-k-1]=='G')
                 sequences[i][start_position[i][j] + k] = 'C';
               if(motifs[j][motif_length-k-1]=='C')
                 sequences[i][start_position[i][j] + k] = 'G';
	       }
	   }
	   else
            for (k = 0; k < motif_length; k++)
            	sequences[i][start_position[i][j] + k] = motifs[j][k];

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
