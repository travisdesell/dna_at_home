cd ../gibbs_cpp/
g++ -O3 -Wall -Wextra motif_models.cpp sampling.cpp sequences.cpp shifting.cpp checkpoint.cpp gibbs_main.cpp ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -lm -o ../bin/gibbs
cd ../bin/
