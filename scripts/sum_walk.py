#!/usr/bin/python
"""identify files in a walk and sum"""
import re
import numpy as np
from optparse import OptionParser
import os


#qsub sum_walk.pbs -v <walk_dir> <walk_num> <min_step>
def main():
    """It all starts here"""
    #usage = "usage: %prog [options] arg"
    usage = "identify files in a walk and sum"
    parser = OptionParser(usage)
    parser.add_option("--min_step", dest="min_step", help="minimum step", metavar="INT", default=None)
    parser.add_option("--walk_num", dest="walk_num", help="walk number", metavar="INT", default=None)
    parser.add_option("--sample_dir", dest="sample_dir", help="What sample should be used?", metavar="STRING", default=None)

    (options, args) = parser.parse_args()
    if not options.min_step or not options.sample_dir or not options.walk_num:
        parser.print_help()
        parser.error("min_step or sample_dir or walk_num not given")

    file_values = get_files_in_walk(options.sample_dir, options.walk_num, int(options.min_step))
    #print file_values
    summed_files = sum_files(file_values)
    print summed_files
    file_out = "%s/walk_%s_steps_-1" % (options.sample_dir, str(options.walk_num))

    np.savetxt(file_out, summed_files, fmt="%s", delimiter=",", newline="\n")

def get_files_in_walk(in_dir, walk_num, min_step):
    """Get files for a walk meeting minimum step"""
    files = sorted(os.listdir(in_dir))

    pattern = re.compile(r"^walk_%s_steps_(\d+)$" % str(walk_num))
    
    file_values = []
    for file_name in files:
        match = pattern.match(file_name)

        if match:
            step = int(match.group(1))
            if step >= min_step:
                file_values.append("%s/%s" % (in_dir, file_name))

    return file_values

def sum_files(file_values):
    """read in each file, sum the results as you go"""

    sum_data = None
    #count = 0
    for in_file in file_values:
        print "summing file: %s" % in_file
        data = None
        with open(in_file, 'r') as in_data:
            data = np.fromstring(in_data.read(), dtype=int, sep=",")

        if sum_data == None:
            sum_data = data
        else:
            sum_data = sum_data + data
        #count += 1
        #if count % 10 == 0:
            #print "count: %s, sum_data: %s" % ((count, sum_data))


    return sum_data


main()
