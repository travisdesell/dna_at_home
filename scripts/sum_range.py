#!/usr/bin/python
"""identify files in a walk and sum"""
import re
from optparse import OptionParser
import os


#qsub sum_walk.pbs -v <walk_dir> <walk_num> <min_step>
def main():
    """It all starts here"""
    #usage = "usage: %prog [options] arg"
    usage = "identify files in a walk and sum"
    parser = OptionParser(usage)
    parser.add_option("--sample_dir", dest="sample_dir", help="What sample should be used?", metavar="STRING", default=None)
    parser.add_option("--no_exec", dest="no_exec", help="Perform the sum?", action="store_true", default=False)

    (options, args) = parser.parse_args()
    if not options.sample_dir:
        parser.print_help()
        parser.error("need sample_dir")

    pattern = re.compile(r"^walk_\d+_steps_range_(\d+)_(\d+)$")
    file_temp = get_files(options.sample_dir, pattern)
    file_values = file_temp[0]
    min_step = file_temp[1]
    total = file_temp[2]
    #print file_values
    summed_files = sum_files(file_values)
    #print summed_files
    file_out = None

    file_out = "%s/walk_steps_all_%s_%s" % (options.sample_dir, str(min_step), str(total))

    print "Output file: %s" % file_out
    #np.savetxt(file_out, summed_files, fmt="%s", delimiter=",", newline="\n")

    with open(file_out, 'w') as out_data:
        for line in summed_files:
            join_string = ",".join([str(val) for val in line])
            print >>out_data, join_string

def get_files(in_dir, pattern):
    """Get files for a walk meeting minimum step"""
    files = sorted(os.listdir(in_dir))
    total = 0
    file_values = []
    min_step = 500000
    for file_name in files:
        match = pattern.match(file_name)

        if match:
            step = int(match.group(1))
            max_step = int(match.group(2))
            size = max_step - step
            if size > 0:
                file_values.append("%s/%s" % (in_dir, file_name))
                if step < min_step:
                    min_step = step
                total += max_step - step

    return [file_values, min_step, total]

def sum_files(file_values):
    """read in each file, sum the results as you go"""

    sum_data = None
    #count = 0
    for in_file in file_values:
        print "summing file: %s" % in_file
        data = []
        with open(in_file, 'r') as in_data:

            for line in in_data:
                line_no_end = line.rstrip('\n')
                list_container = []
                list_container = line_no_end.split(",")
                data.append(list_container)

        if sum_data == None:
            sum_data = data
        else:
            sum_data = [[int(sum_data[y][x]) + int(data[y][x]) for x in range(len(data[y]))] for y in range(len(data))]
        #count += 1
        #if count % 10 == 0:
            #print "count: %s, sum_data: %s" % ((count, sum_data))


    return sum_data

main()
