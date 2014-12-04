#!/usr/bin/python
"""Test for Gibbs Convergence"""
import numpy as np
import time
#from pymc import geweke
from scipy.stats import ks_2samp
import matplotlib as mpl
mpl.use('Agg') #tells pyplot this is a headless node
import matplotlib.pyplot as plt
import os
import re
import sys
import json

USAGE = "\n".join([
    "python ks_intrawalk_convergence.py <IN_DIR> <OUT_GRAPH> <OUT_DATA> [plot_data_file_name]",
    "Plots the kolmogorov Smirnov Statistic which compares two distributions",
    "Each step file within a walk represents a new distribution",
    "Over a number of steps those distributions should start to converge",
    "Displays two figures:",
    "  Distance between CDF's",
    "  A flattening of the intrawalk distance indicates convergence",
    "  Probability that CDF's are from the same distribution",
    "  A high probability indicates convergence",
    "",
    "IN_DIR is the directory containing walks of the form walk_#_steps_#",
    "OUT_GRAPH is the basename of the graph file to output",
    "OUT_DATA is the basename of the data file to output",
    "plot_data_file_name is optional and is used to update a plot with new files",
    "steps that have already been calculated will be loaded"
])

def main():
    """Just so we don't have a global variable mess"""

    if len(sys.argv) < 4:
        print "missing args? : %s, %s" % (sys.argv, USAGE)
        sys.exit(1)

    in_dir = sys.argv[1]
    print "in_dir: %s" % in_dir
    out_graph_file_name = sys.argv[2] + ".png"
    print "out_graph_file_name: %s" % out_graph_file_name
    out_data_file_name = sys.argv[3] + ".json"
    print "out_data_file_name: %s" % out_data_file_name
    plot_data_file_name = None
    if len(sys.argv) == 5:
        plot_data_file_name = sys.argv[4]
        print "plot_data_file_name: %s" % plot_data_file_name

    (file_values, total_files) = get_files_in_dir(in_dir)
    data = process_files_intrawalk(file_values, in_dir, out_data_file_name, total_files, plot_data_file_name)
    plot_data(data, out_graph_file_name)

def ks_2(file_a=None, file_b=None, data_a=None, data_b=None):
    """Perform the Komogorov Smirnov 2 sample comparison"""
    start_time = time.time()
    if data_a is None:
        if file_a is None:
            print "file_a is None and data_a is None"
            raise ValueError()
        with open(file_a, 'r') as thing_a:
            #data_a = np.fromstring(thing_a.read(), dtype=int, sep=",").sort()
            data_a = np.fromstring(thing_a.read(), dtype=int, sep=",")

    if data_b is None:
        if file_b is None:
            raise ValueError("file_b is None and data_b is None")
        with open(file_b, 'r') as thing_b:
            #data_b = np.fromstring(thing_b.read(), dtype=int, sep=",").sort()
            data_b = np.fromstring(thing_b.read(), dtype=int, sep=",")

    data_read_time = time.time()
    dr_time = (data_read_time - start_time)
    ks_data = ks_2samp(data_a, data_b)
    ks_time = time.time() - data_read_time
    return [data_a, data_b, ks_data, dr_time, ks_time]

def get_files_in_dir(in_dir):
    """Organize files by walk and sort"""
    files = sorted(os.listdir(in_dir))

    pattern = re.compile("^walk_(\d+)_steps_(\d+)$")

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


def process_files_intrawalk(file_values, in_dir, out_data_file_name, total_files, plot_data_file_name=None):
    """Given a set of file names, generate KS for them"""
    data = {}
    print "pdfn: %s" % plot_data_file_name
    if plot_data_file_name is not None:
        print "using plot data"
        with open(plot_data_file_name, "r") as plot_data_file:
            data = json.load(plot_data_file)
    #print data

    (out_path, name) = os.path.split(in_dir)
    processed_files = 0
    for walk in file_values:
        walk_str = str(walk)
        #print "data2: %s" % data
        file_name_a = ""
        file_name_b = ""
        data_b = None
        #print "walk: %s" % walk
        #print "walk_str: %s" % walk_str
        if walk_str in data:

            print "found walk in data: %s" % walk_str
        for step in file_values[walk]:
            processed_files += 1
            step_str = str(step)
            print "processing %s of %s" % (processed_files, total_files)
            #print "data3: %s" % data
            if step == 0:
                file_name_b = "_".join(["%s/walk" % in_dir, walk_str, "steps", step_str])
                if walk_str not in data:
                    data[walk_str] = {}
                continue

            #print "data4: %s" % data

            file_name_a = file_name_b
            file_name_b = "_".join(["%s/walk" % in_dir, walk_str, "steps", step_str])
            #print "step: %s" % step
            #print "str step: %s" % str(step)
            #print "step type: %s" % type(step)
            #print "data5: %s" % data
            #print "data[walk]: %s" % data[walk_str]
            #print "keys: %s" % data[walk_str].keys()
            #print "type keys: %s" % type(data[walk_str].keys())
            if step_str in data[walk_str]:
                print "walk: %s, step: %s found data from previous run" % (walk_str, step_str)
                continue

            results = []
            if data_b is not None:
                results = ks_2(file_b=file_name_b, data_a=data_b)
            else:
                results = ks_2(file_a=file_name_a, file_b=file_name_b)
            results.append(file_name_a)
            results.append(file_name_b)

            data_b = results[1]
            data[walk_str][step_str] = results[2:]
             #print results[2:]
    #        print step
    #        if step > 30000:
    #            break
#        if processed_files > 240:#500:
#            break
    with open(out_data_file_name, "w") as out_data_file:
        out_data_file.write(json.dumps(data))
        out_data_file.close()

    return data

def plot_data(data, out_graph_file_name):
    """Generate graphs using matplotlib and pyplot"""
    plt.figure(1)
    #pivot
    for walk in data:
        x1_data = sorted([int(n) for n in data[walk]])
        y1_data = [data[walk][str(step)][0][0] for step in x1_data]
        x2_data = x1_data
        y2_data = [data[walk][str(step)][0][1] for step in x1_data]

        plt.subplot(211)
        plt.ylabel("Distance")
        plt.plot(x1_data, y1_data)
        plt.subplot(212)
        plt.xlabel("step")
        plt.ylabel("Probability")
        plt.plot(x2_data, y2_data)

    plt.savefig(out_graph_file_name)
    #plt.show()
    #print "KS results: %s" % data





if __name__ == "__main__":
    main()

