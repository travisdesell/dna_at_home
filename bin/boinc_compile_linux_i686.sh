cd ../gibbs_cpp/
g++ -O2 -m32 -ftree-vectorize -funroll-loops -Wall -D_BOINC_ -D_DEBUG_ -export-dynamic -I/home/deselt/Software/boinc -I/home/deselt/Software/boinc/api -I/home/deselt/Software/boinc/lib gibbs_main.cpp motif_models.cpp sampling.cpp sequences.cpp shifting.cpp util.cpp phylogeny.cpp checkpoint.cpp ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -o ../bin/Gibbs_$1_i686-pc-linux-gnu -L/home/deselt/Software/boinc/lib -L/home/deselt/Software/boinc/api -lboinc_api -lboinc -pthread
cd ../bin/
