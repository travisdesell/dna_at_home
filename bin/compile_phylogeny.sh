g++ -Wall -DSM_DISALLOW_CACHE -DPHYLOGENY_TEST -I/Users/deselt/Software/gsl-1.15_include /Users/deselt/Software/gsl-1.15_libs/libgslcblas.a /Users/deselt/Software/gsl-1.15_libs/libgsl.a ../gibbs_cpp/phylogeny.cpp ../gibbs_cpp/sequences.cpp ../gibbs_cpp/motif_models.cpp ../../mersenne_twister/dSFMT.c ../gibbs_cpp/lee_phylogeny/SubstitutionMatrix.cc ../gibbs_cpp/lee_phylogeny/SubstitutionModel.cc -o phylogeny

#../gibbs_cpp/lee_phylogeny/Matrix.cc
