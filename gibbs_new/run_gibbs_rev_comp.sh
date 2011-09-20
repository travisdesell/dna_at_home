rm gibbs_sites_checkpoint.txt
rm gibbs_samples_checkpoint.txt
time ./$1 --max_sites 3 --blocks 0.0667 0.4667 0.4667 --motifs 2 reverse_complement,16 reverse_complement,16 --sequence_file $2 --burn_in_period 100000 --sample_period 900000 --enable_shifting 2 5 --print_best_sites 0.5 --print_current_sites --print_accumulated_samples
