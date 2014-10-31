# USAGE EXAMPLE:
# ./print_sampled_sites --motifs forward,6 reverse,6 forward,6 reverse,6 --sequence_file /data/dna_at_home/wgEncodeOpenChromChipMcf7Pol2SerumstimRawDataRep1_deduplicated_100.fa --samples_file /data/dna_at_home/test_hg19_100fa_4/walk_173689_steps_0 --samples_period 10000 --enable_shifting 2 5 --max_sites 4 --best_site_percentage 0.3

g++ ../gibbs_cpp/print_sampled_sites.cpp ../gibbs_cpp/checkpoint.cpp ../gibbs_cpp/motif_models.cpp ../gibbs_cpp/sequences.cpp ../gibbs_cpp/sampling.cpp ../gibbs_cpp/phylogeny.cpp ../gibbs_cpp/util.cpp ../../mersenne_twister/dSFMT.c -DDSFMT_MEXP=19937 -lm -o print_sampled_sites
