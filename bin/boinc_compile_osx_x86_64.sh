cd ../gibbs_cpp/
g++ -fvisibility=hidden -O3 -Wall -D_BOINC_ -DSM_DISALLOW_CACHE -I/Users/deselt/Software/boinc -I/Users/deselt/Software/boinc/api -I/Users/deselt/Software/boinc/lib -L/Users/deselt/Software/boinc/mac_build/build/Deployment -lboinc_api -lboinc gibbs_main.cpp motif_models.cpp sampling.cpp sequences.cpp shifting.cpp util.cpp phylogeny.cpp checkpoint.cpp ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -o ../bin/Gibbs_$1_x86_64-apple-darwin
cd ../bin/
