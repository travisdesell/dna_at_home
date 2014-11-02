import csv
import math
import numpy
import os
import fnmatch      #file name matching
import re           #regular expressions
import string
import sys
from collections import defaultdict


num_motifs = int(sys.argv[1])
directory = sys.argv[2]
walk_id = sys.argv[3]

print "num motifs: ", num_motifs
print "directory: ", directory
print "walk id: ", walk_id


sequence = int(0)
current_motif = int(0)

all_samples = []

for file in sorted(os.listdir(directory)):
    if (fnmatch.fnmatch(file, "walk_" + walk_id + "_*")):
        samples = []

        with open(directory + "/" + file, 'rb') as csvfile:
            samples = [[[] for i in range(num_motifs)]]

            samples_reader = csv.reader(csvfile, delimiter=',', quotechar='#')

            for row in samples_reader:
                if (sequence == len(samples)):
                    samples.append( [[] for i in range(num_motifs)] )

                for sc in row:
                    sample_count = int(sc)
                    #print "samples[{0:4d}][{1:4d}]: {2:4d}".format(sequence, current_motif, sample_count)
                    samples[sequence][current_motif].append(sample_count)

                current_motif = current_motif + 1
                if (current_motif == num_motifs):
                    current_motif = 0
                    sequence = sequence + 1

            #print samples
            print "read file: ", (directory + "/" + file)
        all_samples.append(samples)
        current_motif = 0
        sequence = 0


print "number of files read: ", len(all_samples)

def histogram_difference(samples1,  samples2):
    diff = 0
    min_diff = sys.float_info.max
    max_diff = 0
    count = 0

    for i in range(0, len(samples1)):
        for j in range(0, len(samples1[i])):
            for k in range(0, len(samples1[i][j])):
                c_diff = math.fabs(samples1[i][j][k] - samples2[i][j][k])
                diff = diff + c_diff

                if min_diff > c_diff:
                    min_diff = c_diff
                if max_diff < c_diff:
                    max_diff = c_diff

                count += 1

    print "min_diff: {0}, avg_diff: {1} max_diff: {2}".format(min_diff, diff / count, max_diff),

    return diff

for correlation in range(1, 11):
    first = 0
    avg_distance = 0.0
    print "correlation {0}:".format(correlation)
    while (first + correlation < len(all_samples)):
        distance = histogram_difference(all_samples[first], all_samples[first+correlation])
        avg_distance += distance
        print " - total_distance: {0}".format(distance)
        first = first + 1

    avg_distance /= len(all_samples) - correlation
    print "\naverage: {0}\n".format(avg_distance)

