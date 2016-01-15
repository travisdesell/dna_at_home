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
walk_depth = sys.argv[3]

print "num motifs: ", num_motifs
print "directory: ", directory
print "walk_depth: ", walk_depth


sequence = int(0)
current_motif = int(0)

all_samples = []

for file in sorted(os.listdir(directory)):
    if (fnmatch.fnmatch(file, "walk_*_" + walk_depth)):
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
    for i in range(0, len(samples1)):
        for j in range(0, len(samples1[i])):
            for k in range(0, len(samples1[i][j])):
                diff = diff + math.fabs(samples1[i][j][k] - samples2[i][j][k])
    return diff

distances = numpy.zeros((len(all_samples), len(all_samples)), dtype=numpy.float64)
print distances

for walk1 in range(0, len(all_samples)):
    for walk2 in range (walk1 + 1, len(all_samples)):
        distances[walk1][walk2] = histogram_difference(all_samples[walk1], all_samples[walk2])
        distances[walk2][walk1] = distances[walk1][walk2]

    min = sys.float_info.max
    max = 0
    mean = 0
    for i in range(0, len(all_samples)):
        if i == walk1:
            continue
        
        if distances[walk1][i] < min:
            min = distances[walk1][i]
        if distances[walk1][i] > max:
            max = distances[walk1][i]
        mean += distances[walk1][i];


    mean /= len(all_samples) - 1
    print "walk {0}, min: {1}, mean {2}: max {3}".format( walk1, min, mean, max )
    overall_mean += mean
print distances

overall_mean /= len(all_samples)

print"\noverall mean: {0}\n".format(overall_mean)

