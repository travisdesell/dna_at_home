rm gibbs_sites_checkpoint.txt
rm gibbs_samples_checkpoint.txt
rm accumulated_samples.txt
time ./$1 --max_sites 3 --blocks 0.1 0.3 0.3 0.3 --motifs palindromic,22 palindromic,22 --sequence_file $2 --burn_in_period 100000 --sample_period 600000 --enable_shifting 2 5 --print_best_sites 0.005 --print_accumulated_samples --phylogeny_file $3 --tt_factor 3.0
