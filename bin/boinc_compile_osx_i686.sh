cd ../gibbs_new/
/usr/bin/g++-4.0 -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -mmacosx-version-min=10.4 -O3 -msse3 -Wall -mfpmath=sse -mtune=prescott -ftree-vectorize -funroll-loops -D_BOINC_ -I/Users/deselt/Software/boinc -I/Users/deselt/Software/boinc/api -I/Users/deselt/Software/boinc/lib -L/Users/deselt/Software/boinc/mac_build/build/Deployment -lboinc_api -lboinc gibbs_main.c motif_models.c sampling.c sequences.c shifting.c checkpoint.c ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -o ../bin/Gibbs_$1_i686-apple-darwin
cd ../bin/
