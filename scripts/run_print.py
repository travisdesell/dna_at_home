#!/usr/bin/python
"""create print_sampled_sites_command"""
import MySQLdb
import re
from optparse import OptionParser
import subprocess
import shlex
import os
import sys

def main():
    """where it all begins"""
    usage = "Generate new print_sampled_sites commands which in turn generate .motif files from walk files"
    parser = OptionParser(usage)
    #parser.add_option("--step", dest="step", help="Which step should be compared?", metavar="INT")
    parser.add_option("--sample_name", dest="sample_name", help="What sample should be used?", metavar="STRING")
    parser.add_option("--passwd", dest="passwd", help="What passwd to use", metavar="STRING")
    parser.add_option("--best_pct", dest="best_pct", help="What is the minimum representation required", metavar="STRING")
    parser.add_option("--dna_root", dest="dna_root", help="where is the dna_at_home directory?", metavar="STRING")
    parser.add_option("--data_root", dest="data_root", help="where is the data directory?", metavar="STRING")
    parser.add_option("--seq_root", dest="seq_root", help="where is the seq directory?", metavar="STRING")
    parser.add_option("--db_host", dest="db_host", help="db_host", metavar="STRING")
    parser.add_option("--no_exec", dest="no_exec", help="writeout a shell script", action="store_true", default=False)

    (options, args) = parser.parse_args()
    if not options.sample_name:
        parser.print_help()
        parser.error("step or sample_name not given")

    #print "options: %s" % options

    database = MySQLdb.connect(
        host=options.db_host,
        user="tdesell",
        passwd=options.passwd,
        port=3306,
        db="csg"
    )

    cur = database.cursor()
    #print """SELECT name, sequences_filename, command_line_options FROM gibbs_sampler WHERE name like %s""" % (options.sample_name,)
    cur.execute(
        """SELECT name, sequences_filename, command_line_options FROM gibbs_sampler WHERE name like %s""",
        (options.sample_name,)
    )
    
    result = cur.fetchone()
    #print "result: %s" % (result,)
    fields = map(lambda x: x[0], cur.description)
    #result = [dict(zip(fields, row)) for row in cur.fetchall()]
    sampler = dict(zip(fields, result))
    #print "sampler: %s" % (sampler,)
#    if result == 0:
#        parser.print_help()
#        parser.error("sample_name not found")


    # --max_sites 3 --blocks 0.1 0.3 0.3 0.3 --motifs  forward,6 reverse,6 --enable_shifting 2 5 --print_best_sites 0.1 --checkpoint_frequency 10000 --sequence_file sequences.txt --burn_in_period 0 --sample_period 10000
    max_sites_r = re.compile(r"--max_sites \d+")
    shift_r = re.compile(r"--enable_shifting( \d+)+")
    best_r = re.compile(r"--print_best_sites (\d*\.\d+)")
    #sample_period_r = re.compile(r"--sample_period (\d+)")
    #NOTE samples_period has been changed since we combined all data into one file
    #use the last number in the file name
    #sample_period_r = re.compile(r"--sample_period (\d+)")


    #./print_sampled_sites --max_sites 4 --motifs forward,6 reverse,6 forward,6 reverse,6  --enable_shifting 2 5 --best_site_percentage .1 --sequence_file /data/dna_at_home/snail/wgEncodeOpenChromChipMcf7Pol2SerumstimRawDataRep1_deduplicated_snail_1000.fa --samples_period 10000    --samples_file /data/dna_at_home/snail_hg19_1000fa_1/walk_181178_steps_250000 > /data/dna_at_home/snail_hg19_1000fa_1/walk_181178_steps_250000.motifs
    #sampler: {'command_line_options': ' --max_sites 4 --blocks  0.1 0.225 0.225 0.225 0.225 --motifs forward,6 reverse,6 forward,6 reverse,6  --enable_shifting 2 5 --print_best_sites 0.1 --checkpoint_frequency 1006 --sequence_file sequences.txt --burn_in_period 0 --sample_period 10000', 'name': 'snail_hg19_1000fa_1'}

    #motif_r = re.compile("--motifs\s*(\s\+\w+,\d+)+")
    motif_r = re.compile(r"--motifs(\s+\w+,\d+)+")

    motifs_raw = motif_r.search(sampler["command_line_options"]).group(0)
    motifs = "--motifs '%s'" % motifs_raw.replace("--motifs ","")
    max_sites = max_sites_r.search(sampler["command_line_options"]).group(0)
    shift = shift_r.search(sampler["command_line_options"]).group(0)
    best = None
    if options.best_pct is None or options.best_pct == "":
        best_m = best_r.search(sampler["command_line_options"])
        best = "--best_site_percentage %s" % best_m.group(1)
    else:
        best = options.best_pct

    #sample_period_m = sample_period_r.search(sampler["command_line_options"])
    #NOTE replaced with value from filename
    #samples_period = "--samples_period %s" % sample_period_m.group(1)
    samples_period = -1

    pattern = re.compile(r"^walk_steps_all_(\d+)_(\d+)$")
    samples_files = get_files("%s/%s" % (options.data_root, sampler["name"]), pattern)
    samples_file = ""
    if len(samples_files[0]) == 0:
        print "no walk_steps_all found"
        sys.exit(0)
    elif len(samples_files[0]) == 1:
        samples_file = samples_files[0][0]
        samples_period = "--samples_period %s" % samples_files[1][0]
    elif len(samples_files[0]) > 1:
        print "found multiple candidates: %s" % samples_files
        sys.exit(0)

    #amples_file = "%s/%s/walk_steps_all_%s" % (options.data_root, sampler["name"], options.step)
    sequence = "--sequence_file %s/%s" % (options.seq_root, os.path.basename(sampler["sequences_filename"]))

    #args = " ".join(("/home/tdesell/dna_at_home/bin/print_sampled_sites", max_sites, motifs, shift, best, samples_period, sequence, "--samples_file %s" % samples_file, ">","%s_motifs/%s.motifs" % (head, tail)))
    args = " ".join(("%s/bin/print_sampled_sites" % options.dna_root, max_sites, motifs, shift, best, samples_period, sequence, "--samples_file %s" % samples_file))

    new_args = shlex.split(args)
    if options.no_exec:
        print "qsub %s/scripts/run.pbs -v SCRIPT=\"%s > %s.motifs\"" % (options.dna_root, args, samples_file)
    else:
        with open("%s.motifs" % samples_file, "w+") as out_file:
            subprocess.call(new_args, stdout=out_file)

def get_files(in_dir, pattern):
    """Get files for a walk meeting minimum step"""
    files = sorted(os.listdir(in_dir))
    #print files
    file_values = []
    periods = []
    for file_name in files:
        match = pattern.match(file_name)

        if match:
            file_values.append("%s/%s" % (in_dir, file_name))
            periods.append(match.group(2))

    return [file_values, periods]


main()
