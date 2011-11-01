#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "../../mersenne_twister/dSFMT.h"

#define MODEL_TYPE_NORMAL 0
#define MODEL_TYPE_REVERSE 1 //Reverse Complement
//#define MODEL_TYPE_FORWARD 2 //Not implemented
#define MODEL_TYPE_PALINDROMIC 3

#define MODEL_BUILD_PROBABILISTIC 0
#define MODEL_BUILD_EXACT 1

#define ALPHABET_LENGTH 4

#define A_POSITION 0
#define C_POSITION 1
#define G_POSITION 2
#define T_POSITION 3

//matrix[motif_width][ALPHABET_LENGTH]
float gaussRandomNum()
{
//Ref:: 'Numerical Recipes in C'
//Returns a normally distributed deviate with zero mean and unit variance
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;
	if (iset == 0) { 
	do {
	v1=2.0*(float)dsfmt_gv_genrand_close_open()-1.0; 
	v2=2.0*(float)dsfmt_gv_genrand_close_open()-1.0; 
	rsq=v1*v1+v2*v2; 
	} while (rsq >= 1.0 || rsq == 0.0);
	fac=sqrt(-2.0*log(rsq)/rsq);
	gset=v1*fac;
	iset=1;
	return v2*fac;
	} else {
	iset=0;
	return gset;
}

}
void fillProbabilities(float **matrix_,int motif_lenght_)
{
	int i;
	float xA,xC,xG,xT,z;
	xA=0;
	xC=0;
	xG=0;
	xT=0;
	z =0;
	//printf("Calculating Probabilities %d\n",motif_lenght_);
	for(i=0;i<motif_lenght_;i++)
	{
	xA=pow(gaussRandomNum(),2);
	xC=pow(gaussRandomNum(),2);
	xG=pow(gaussRandomNum(),2);
	xT=pow(gaussRandomNum(),2);
	z=xA+xC+xG+xT;	
	
	matrix_[i][A_POSITION]=xA/z;
	matrix_[i][C_POSITION]=xC/z;
	matrix_[i][G_POSITION]=xG/z;
	matrix_[i][T_POSITION]=xT/z;
	//printf("END Calculating Probabilities %d\n",i);
	}
	//printf("END Calculating Probabilities \n");
}

int main(int number_arguments, char **arguments) {
    int character, i, j, k;
    char **sequences, **motifs;
    int **start_position;
    float **probabilities, randNum;
    int number_sequences, sequence_length,palindrome_length, insert_all_motifs;
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
    randNum=0;
    insert_all_motifs= 1;
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

        } else if ( !strcmp(arguments[i], "--insert_all_motifs")) {
            insert_all_motifs = atoi(arguments[++i]);

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

    probabilities = (float**)malloc((motif_length) * sizeof(float*));
    for (i = 0; i < motif_length; i++) {
        probabilities[i] = (float*)malloc(ALPHABET_LENGTH * sizeof(float));
        for (j = 0; j < ALPHABET_LENGTH; j++) {
            probabilities[i][j] = 0.0;
        }
    }

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

    if(MODEL_BUILD_PROBABILISTIC==motif_build)
      {
/////////////////////////////////////// START-PROBABILISTIC/////////////////////
	fillProbabilities(probabilities,motif_length);
	printf("*****Probabilities***** \n");
	printf("A	,C	,G	,T \n");
	for (i = 0; i < motif_length; i++) {
        for (j = 0; j < ALPHABET_LENGTH; j++) {
             printf("%f,",probabilities[i][j]);
        }
    	printf("\n");			   }
	
	if(MODEL_TYPE_REVERSE == motif_model||MODEL_TYPE_NORMAL == motif_model)
     {
    	for (i = 0; i < number_motifs; i++) {
        for (j = 0; j < motif_length; j++) {
                randNum = (float)dsfmt_gv_genrand_open_close();

	float uA=probabilities[j][A_POSITION];
        float uG=probabilities[j][A_POSITION]+probabilities[j][G_POSITION];
        float uC=probabilities[j][A_POSITION]+probabilities[j][G_POSITION]+probabilities[j][C_POSITION];

        /*	printf("%f LIMITS-A %f \n",0.0,uA);
		printf("%f LIMITS-G %f \n",uA,uG);
		printf("%f LIMITS-C %f \n",uG,uC);
		printf("%f LIMITS-T %f \n",uC,1.0);
	        printf("j=%d randNum=%f\n",j,randNum);*/
            if (0.0 <= randNum && randNum< uA) motifs[i][j] = 'A';
            if (uA <= randNum && randNum < uG) motifs[i][j] = 'G';
            if (uG <= randNum && randNum < uC) motifs[i][j] = 'C';
            if (uC <= randNum && randNum < 1.0) motifs[i][j] = 'T';
        }
        motifs[i][motif_length] = '\0';
        printf("motif [%d]: [%s]\n", i, motifs[i]);
  	  }
     }

	if(MODEL_TYPE_PALINDROMIC == motif_model)
     {
	if( 2*palindrome_length > motif_length)
             {
                printf("Palindrome lenght should be atmost half of the motif_length\n Palindrome Length=%d Motif Length=%d\n",palindrome_length, motif_length);
                return;
             }

    	for (i = 0; i < number_motifs; i++) {
        for (j = 0; j < motif_length-palindrome_length; j++) {
                randNum = (float)dsfmt_gv_genrand_open_close();

	float uA=probabilities[j][A_POSITION];
        float uG=probabilities[j][A_POSITION]+probabilities[j][G_POSITION];
        float uC=probabilities[j][A_POSITION]+probabilities[j][G_POSITION]+probabilities[j][C_POSITION];

        /*	printf("%f LIMITS-A %f \n",0.0,uA);
		printf("%f LIMITS-G %f \n",uA,uG);
		printf("%f LIMITS-C %f \n",uG,uC);
		printf("%f LIMITS-T %f \n",uC,1.0);
	        printf("j=%d randNum=%f\n",j,randNum);*/
	
	   if (0.0 <= randNum && randNum< uA)
            	{
            	motifs[i][j] = 'A';
		if(j<palindrome_length)
            	motifs[i][motif_length-j-1] = 'T';
            	}
            if (uA <= randNum && randNum < uG)
		{
            	motifs[i][j] = 'G';
		if(j<palindrome_length)
               	motifs[i][motif_length-j-1] = 'C';
		}
            if (uG <= randNum && randNum < uC)
		{
            	motifs[i][j] = 'C';
		if(j<palindrome_length)
                motifs[i][motif_length-j-1] = 'G';
            	}

            if (uC <= randNum && randNum < 1.0)
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
/////////////////////////////////////// END-PROBABILISTIC///////////////////////	
      }

    if(MODEL_BUILD_EXACT==motif_build)
{
///////////////////////////////////START-EXACT//////////////////////////////////
    if(MODEL_TYPE_REVERSE == motif_model||MODEL_TYPE_NORMAL == motif_model)
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
    if(MODEL_TYPE_PALINDROMIC == motif_model)
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
///////////////////////////////////END-EXACT////////////////////////////////////
}
//    printf("overrepresented sequence: [%s]\n", motifs);

    for (i = 0; i < number_sequences; i++) {
	 if(insert_all_motifs ==0)
	 character = (int)(dsfmt_gv_genrand_open_close() * (number_motifs))+1; //select randomly the number of motifs to insert in each sequence
	else character = number_motifs;
	 
        for (j = 0; j < character; j++) {
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
	    if(start_position[i][j]>0)	
            printf(", start pos [%d] subseq [%d][%s] end pos [%d]", start_position[i][j], j, motifs[j], start_position[i][j] + motif_length - 1);
        }
        printf("\n");

        printf("%s\n\n", sequences[i]);
    }

}
