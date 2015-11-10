#!/usr/bin/python
"""Identify walks to sum"""
#import MySQLdb
import re
from optparse import OptionParser
#import subprocess
#import shlex
import os


def main():
    """It all starts here"""
    #usage = "usage: %prog [options] arg"
    usage = "identify walks to sum, print or qsub sum commands"
    parser = OptionParser(usage)
    parser.add_option("--min_step", dest="min_step", help="minimum step", metavar="INT", default=None)
    parser.add_option("--sample_dir", dest="sample_dir", help="What sample should be used?", metavar="STRING", default=None)
    #parser.add_option("--passwd", dest="passwd", help="What passwd to use", metavar="STRING")
    parser.add_option("--no_exec", dest="no_exec", help="writeout a shell script", action="store_true", default=False)

    (options, args) = parser.parse_args()
    if not options.min_step or not options.sample_dir:
        parser.print_help()
        parser.error("min_step or sample_dir not given")

    (file_values, total_files) = get_files_in_dir(options.sample_dir)
    commands = make_commands(file_values, options.sample_dir, options.min_step)
    print commands

    #if options.no_exec:
    #    print "%s > %s.motifs" % (args, samples_file)
    #else:
    #    #with open("%s.motifs" % samples_file, "w+") as out_file:
    #    #    subprocess.call(new_args, stdout=out_file)
    #    pass

def get_files_in_dir(in_dir):
    """Organize files by walk and sort"""
    files = sorted(os.listdir(in_dir))

    pattern = re.compile(r"^walk_(\d+)_steps_(\d+)$")

    file_values = {}
    total_files = 0
    for file_name in files:
        #print file_name
        match = pattern.match(file_name)
        #print pattern
        #print match
        if match:
            total_files += 1
            match1 = int(match.group(1))
            match2 = int(match.group(2))
            if match1 not in file_values:
                file_values[match1] = []
            file_values[match1].append(match2)
    #print file_values

    for file_name in file_values:
        file_values[file_name] = sorted(file_values[file_name])

    #print file_values
    return (file_values, total_files)

def make_commands(file_values, sample_dir, min_step):
    """build some commands"""
    commands = []
    for walk_num in file_values.keys():
        print "WALK_NUM: %s\n" % walk_num

        #qsub sum_walk.pbs -v <walk_dir> <walk_num> <min_step>
        commands.append(" ".join(("qsub", "sum_walk.pbs", "-v", "ARGS=\"%s\"" % " ".join(( sample_dir, str(walk_num), str(min_step))))))
	    #new_args = shlex.split(args) do we need shlex?
    return commands


main()
