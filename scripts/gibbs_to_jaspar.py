#`/usr/bin/python
"""gibbs_to_jaspar converts gibbs motif pwm output to jaspar format by removing all the cruft"""

import sys
import os
import re
from optparse import OptionParser


def process_file(in_name):
    """process a single file"""

    solutions = []
    count_line = re.compile(r"^letter count \[([ACGT])\]:\s+(.*)$")
    body_count = 0
    with open(in_name, 'r') as data:
        header = ">" + in_name
        body = ""

        for line in data:
            line = line.rstrip()
            if line == "<motif_models>":
                continue
            elif line.startswith('Model type: '):
                header += "\t" + line[12:]
            elif line.startswith('Model length: '):
                header += "\t" + line[14:] + "\n"
            else:
                match_count = count_line.match(line)
                if match_count:
                    body += "%s [%s]\n" % (match_count.group(1), match_count.group(2))
                    body_count += 1
                    if body_count == 4:
                        solutions.append("%s%s" % (header, body))
                        body_count = 0
                        header = ">" + in_name
                        body = ""
                else:
                    continue

        #for solution in solutions:
        #    print solution
        return solutions

def get_files_in_dir(in_dir):
    """Organize files by walk and sort.  borrowed from ks_interwalk_convergence.py searching for motifs isntead of samples here"""

    parents = sorted(os.listdir(in_dir))
    print parents
    file_values = {}
    for parent in parents:
        files = []
        try:
            files = sorted(os.listdir("%s/%s" % (in_dir, parent)))
        except Exception as ex:
            #print ex
            continue

        pattern = re.compile(r"^motifs_(\d+)_(\d+)$")

        total_files = 0
        for file_name in files:
            #print file_name
            match = pattern.match(file_name)
            #print pattern
            if match:
                total_files += 1
                match1 = int(match.group(1))
                match2 = int(match.group(2))
                if match1 not in file_values:
                    file_values[match1] = []
                file_values[match1].append(match2)
            else:
                print "nomatch:", file_name
    #print file_values

    for walk_name in file_values:
        file_values[walk_name] = sorted(file_values[walk_name])

    #print file_values
    return (file_values, total_files)



def main():
    """keep it all organized"""
    usage = "USAGE: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--in_dir", dest="in_dir", help="a directory containing walk directories", default=None, metavar="STRING")

    accumulations = {}
    (options, args) = parser.parse_args()
    #print "OPTIONS: ", options
    #print "ARGS: ", args


    (file_values, total_files) = get_files_in_dir(options.in_dir)

    print "file_values:", file_values
    for walk in file_values:
        for step in file_values[walk]:
            if step not in accumulations:
                accumulations[step] = {"forward": [], "reverse": []}

            file_data = process_file("walk_%s/motifs_%s_%s" % (walk, walk, step))
            for index in range(len(file_data)):
                if index % 2 == 0:
                    #even, forward
                    accumulations[step]["forward"].append(file_data[index])
                else:
                    accumulations[step]["reverse"].append(file_data[index])

    for step in accumulations:
        for direction in accumulations[step]:
            (head, tail) = os.path.split(os.path.realpath(options.in_dir))
            with open("%s/%s/motifs_%s_%s_%s.txt" % (head, tail, tail, step, direction), "w+") as out_file:
                for motif in accumulations[step][direction]:
                    out_file.write(motif)
    #exit()

if __name__ == "__main__":
    main()
