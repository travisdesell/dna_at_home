#!/usr/bin/python
import MySQLdb
import re
from optparse import OptionParser
import subprocess
import shlex
import os


#usage = "usage: %prog [options] arg"
usage = ""
parser = OptionParser(usage)
parser.add_option("--step", dest="step", help="Which step should be compared?", metavar="INT")
parser.add_option("--sample_name", dest="sample_name", help="What sample should be used?", metavar="STRING")
parser.add_option("--passwd", dest="passwd", help="What passwd to use", metavar="STRING")
parser.add_option("--best_pct", dest="best_pct", help="Whatis the minimum represnetation required", metavar="STRING")
parser.add_option("--dna_root", dest="dna_root", help="where is the dna_at_home directory?", metavar="STRING")
parser.add_option("--data_root", dest="data_root", help="where is the data directory?", metavar="STRING")
parser.add_option("--seq_root", dest="seq_root", help="where is the seq directory?", metavar="STRING")
parser.add_option("--no_exec", dest="no_exec", help="writeout a shell script", action="store_true", default=False)

(options, args) = parser.parse_args()
if not options.step or not options.sample_name:
	parser.print_help()
	parser.error("step or sample_name not given")


db = MySQLdb.connect(
		host="localhost", 
		user="tdesell",
		passwd=options.passwd,
		db="wildlife")

cur = db.cursor()
cur.execute(
	"\n".join((
		"SELECT name, sequences_filename, command_line_options, gibbs_walk.id, gibbs_walk.current_steps",
		"FROM gibbs_sampler JOIN gibbs_walk ON sampler_id = gibbs_sampler.id",
		"WHERE name like %s"
		)),
	options.sample_name
)

fields = map(lambda x:x[0], cur.description)
result = [dict(zip(fields, row)) for row in cur.fetchall()]
if len(result) == 0:
	parser.print_help()
	parser.error("sample_name not found")

#motif_r = re.compile("--motifs\s*(\s\+\w+,\d+)+")
motif_r = re.compile("--motifs(\s+\w+,\d+)+")

# --max_sites 3 --blocks 0.1 0.3 0.3 0.3 --motifs  forward,6 reverse,6 --enable_shifting 2 5 --print_best_sites 0.1 --checkpoint_frequency 10000 --sequence_file sequences.txt --burn_in_period 0 --sample_period 10000
max_sites_r = re.compile("--max_sites \d+")
shift_r = re.compile("--enable_shifting( \d+)+")
best_r = re.compile("--print_best_sites (\d*\.\d+)")
sample_period_r = re.compile("--sample_period (\d+)")
#_r = re.compile("--print_best_sites (\d*\.\d+)")


#./print_sampled_sites --max_sites 4 --motifs forward,6 reverse,6 forward,6 reverse,6  --enable_shifting 2 5 --best_site_percentage .1 --sequence_file /data/dna_at_home/snail/wgEncodeOpenChromChipMcf7Pol2SerumstimRawDataRep1_deduplicated_snail_1000.fa --samples_period 10000    --samples_file /data/dna_at_home/snail_hg19_1000fa_1/walk_181178_steps_250000 > /data/dna_at_home/snail_hg19_1000fa_1/walk_181178_steps_250000.motifs
#WALK: {'current_steps': 500000L, 'command_line_options': ' --max_sites 4 --blocks  0.1 0.225 0.225 0.225 0.225 --motifs forward,6 reverse,6 forward,6 reverse,6  --enable_shifting 2 5 --print_best_sites 0.1 --checkpoint_frequency 1006 --sequence_file sequences.txt --burn_in_period 0 --sample_period 10000', 'name': 'snail_hg19_1000fa_1', 'id': 180511L}

for walk in result:
#	print "WALK: %s\n" % walk

#	print "WALK_CMD: %s\n" % walk["command_line_options"]
	motifs = motif_r.search(walk["command_line_options"]).group(0)
	#print "motifs: '%s'" % motifs
	max_sites = max_sites_r.search(walk["command_line_options"]).group(0)
	#print "max_sites: '%s'" % max_sites
	shift = shift_r.search(walk["command_line_options"]).group(0)
	best = None
	if options.best_pct is None or options.best_pct == "":
		best_m = best_r.search(walk["command_line_options"])
		best = "--best_site_percentage %s" % best_m.group(1)
	else:
		best = options.best_pct
	#print "best: '%s'" % best
	sample_period_m = sample_period_r.search(walk["command_line_options"])
	samples_period = "--samples_period %s" % sample_period_m.group(1)
	#print "samples_period: '%s'" % samples_period
	
	samples_file = "%s/%s/walk_%s_steps_%s" % (options.data_root, walk["name"], walk["id"], options.step)
	#print "samples_file: '%s'" % samples_file

	sequence = "--sequence_file %s" % walk["sequences_filename"]
	sequence = "%s/%s" % (options.seq_root, os.path.basename(sequence))


	#XXX fix the hardcoded path
	#args = " ".join(("/home/tdesell/dna_at_home/bin/print_sampled_sites", max_sites, motifs, shift, best, samples_period, sequence, "--samples_file %s" % samples_file, ">","%s_motifs/%s.motifs" % (head, tail)))
	args = " ".join(("%s/bin/print_sampled_sites" % options.dna_root, max_sites, motifs, shift, best, samples_period, sequence, "--samples_file %s" % samples_file))
#	print "args: %s\n" % args
#	print "real: /home/tdesell/dna_at_home/bin/print_sampled_sites --max_sites 4 --motifs forward,6 reverse,6 forward,6 reverse,6 --enable_shifting 2 5 --best_site_percentage 0.1 --samples_period 10000 --sequence_file /data/dna_at_home/snail/wgEncodeOpenChromChipMcf7Pol2SerumstimRawDataRep1_deduplicated_snail_1000.fa --samples_file /data/dna_at_home/snail_hg19_1000fa_1/walk_181178_steps_250000 > /data/dna_at_home/snail_hg19_1000fa_1/walk_181178_steps_250000.motifs\n"
	new_args = shlex.split(args)
#	print "new : %s\n" % new_args
	if options.no_exec:
		print "%s\n" % args
	else:
		with open("%s.motifs" % samples_file, "w+") as out_file: 
			subprocess.call(new_args, stdout=out_file)
		

#	steps = 
	



