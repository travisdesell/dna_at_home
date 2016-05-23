#rm gibbs_sites_checkpoint.txt
#rm gibbs_samples_checkpoint.txt
time ./$1 --max_sites 2 --blocks 0.0667 0.4667 0.4667 --motifs forward,16 reverse,16 --sequence_file $2 --burn_in_period 0 --sample_period 1000 --enable_shifting 2 5 --print_best_sites 0.1 --seed $3 --checkpoint_frequency 500 --current_sites input_sites.txt --print_motif_models
#time ./$1 --max_sites 2 --blocks 0.0667 0.4667 0.4667 --motifs forward,16 reverse,16 --sequence_file $2 --burn_in_period 0 --sample_period 5000 --enable_shifting 2 5 --print_best_sites 0.1 --print_current_sites --print_accumulated_samples --seed 12345 --checkpoint_frequency 1000
