rm gibbs_sites_checkpoint.txt
rm gibbs_samples_checkpoint.txt
time ./$1 --max_sites 2 --blocks 0.0667 0.4667 0.4667 --motifs 2 forward,16 reverse,16 --sequence_file $2 --burn_in_period 10000 --sample_period 0 --enable_shifting 2 5 --print_best_sites 0.1 --print_current_sites --print_accumulated_samples
