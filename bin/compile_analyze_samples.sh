cd ../gibbs_cpp/
#g++ -g -rdynamic -ansi -D_DEBUG_ -O3 -Wall -Wextra -fno-nonansi-builtins -pedantic-errors -Wall -Wextra -Wshadow -Wno-unused-parameter -Wcast-qual -Wcast-align -Wredundant-decls -Wabi -Wctor-dtor-privacy -Weffc++ -Wold-style-cast -Woverloaded-virtual -Werror motif_models.cpp sampling.cpp sequences.cpp shifting.cpp checkpoint.cpp gibbs_main.cpp ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -lm -o ../bin/gibbs_cpp
g++ -g -rdynamic -D_DEBUG_ -O3 -Wall -Wextra motif_models.cpp sampling.cpp sequences.cpp shifting.cpp checkpoint.cpp analyze_samples.cpp ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -lm -o ../bin/analyze_samples
#g++ -O3 -Wall -Wextra motif_models.cpp sampling.cpp sequences.cpp shifting.cpp checkpoint.cpp gibbs_main.cpp ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -lm -o ../bin/gibbs_cpp
cd ../bin/
