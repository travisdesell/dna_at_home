cd ../gibbs_new/
g++ -O3 -Wall -D_BOINC_ -I/Users/deselt/Software/boinc -I/Users/deselt/Software/boinc/api -I/Users/deselt/Software/boinc/lib -L/Users/deselt/Software/boinc/mac_build/build/Deployment -lboinc_api -lboinc gibbs_main.c motif_models.c sampling.c sequences.c shifting.c checkpoint.c ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -o ../bin/Gibbs_$1_x86_64-apple-darwin
cd ../bin/
