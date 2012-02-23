cd ../gibbs_cpp/
g++ -fvisibility=hidden -O3 -Wall -D_BOINC_ -DSM_DISALLOW_CACHE -I/home/boinc/Software/boinc -I/home/boinc/Software/boinc/api -I/home/boinc/Software/boinc/lib gibbs_main.cpp motif_models.cpp sampling.cpp sequences.cpp shifting.cpp util.cpp phylogeny.cpp checkpoint.cpp ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -o ../bin/Gibbs_$1_x86_64-pc-linux-gnu -L/home/boinc/Software/boinc/lib -L/home/boinc/Software/boinc/api -lboinc_api -lboinc -pthread
cd ../bin/
