cd ../gibbs_new/
g++ -O3 -funroll-loops -Wall -D_BOINC_ -I/home/deselt/software/boinc -I/home/deselt/software/boinc/api -I/home/deselt/software/boinc/lib gibbs_main.c motif_models.c sampling.c sequences.c shifting.c checkpoint.c ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -o ../bin/Gibbs_$1_x86_64-pc-linux-gnu -L/home/deselt/software/boinc/lib -L/home/deselt/software/boinc/api -lboinc_api -lboinc -pthread
cd ../bin/
