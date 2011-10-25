cd ../gibbs_new/
gcc -O3 -Wall -Wextra gibbs_main.c motif_models.c sampling.c sequences.c shifting.c checkpoint.c ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -lm -o ../bin/gibbs
cd ../bin/
