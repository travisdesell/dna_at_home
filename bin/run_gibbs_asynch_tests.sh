rm gibbs_sites_checkpoint.txt
rm gibbs_samples_checkpoint.txt
time ./$1 --max_sites 2 --blocks 0.0667 0.4667 0.4667 --motifs forward,16 reverse,16 --sequence_file $2 --burn_in_period 16000 --print_current_sites_frequency 500 --print_sites_frequency_logarithmic --total_independent_walks 3 --seed 5 --enable_shifting 2 5
