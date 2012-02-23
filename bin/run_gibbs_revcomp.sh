#rm accumulated_samples.txt
#rm gibbs_sites_checkpoint.txt
#rm gibbs_samples_checkpoint.txt
time ./$1 --max_sites 4 --blocks 0.1 0.225 0.225 0.225 0.225 --motifs forward,16 reverse,16 forward,16 reverse,16 --sequence_file $2 --burn_in_period 100000 --sample_period 500000 --print_best_sites 0.005 --print_accumulated_samples --enable_shifting 2 5
