#!/usr/bin/python
"""Test for Gibbs Convergence"""

import sys
import os
#print sys.path
#print os.listdir("/usr/lib64/python2.6/site-packages")
import numpy as np
import time
import math
#from pymc import geweke
from scipy.stats import ks_2samp
import matplotlib as mpl
mpl.use('Agg') #tells pyplot this is a headless node
import matplotlib.pyplot as plt
import re
import json
import MySQLdb




USAGE = "\n".join([
    "python ks_intrawalk_convergence.py <IN_DIR> <OUT_GRAPH> db_passwd",
    #"$ ./ks_convergence_intrawalk.py /home/kzarns/test_hg19_10fa_5 test_graph6 passwd",
    "$ ./ks_convergence_intrawalk.py /home/kzarns/test_hg19_10fa_5 \"graph title\" passwd",
    "Plots the kolmogorov Smirnov Statistic which compares two distributions",
    "Each step file within a walk represents a new distribution",
    "Over a number of steps those distributions should start to converge",
    "Displays two figures:",
    "  Distance between CDF's",
    "  A flattening of the intrawalk distance indicates convergence",
    "  Probability that CDF's are from the same distribution",
    "  A high probability indicates convergence",
    "",
    "  Distance and probability stats",
    "IN_DIR is the directory containing walks of the form walk_#_steps_#",

    #    "OUT_GRAPH is the basename of the graph file to output",
    #    "OUT_DATA is the basename of the data file to output",
    #    "plot_data_file_name is optional and is used to update a plot with new files",
    #    "steps that have already been calculated will be loaded"
])

def main():
    """Just so we don't have a global variable mess"""




    in_dir = None
    out_graph_file_name = None
    out_data_file_name = None
    plot_data_file_name = None
    if len(sys.argv) != 4:
        print "bad dir? : %s, %s" % (sys.argv, USAGE)
        sys.exit(1)

    in_dir = sys.argv[1]
    print "in_dir: %s" % in_dir

    title = sys.argv[2]

    passwd = sys.argv[3]
    (head, tail) = os.path.split(in_dir)
    if tail == "":
        (head, tail) = os.path.split(head)

    out_graph_file_name = "./%s_full.png" % (tail)
    out_stats_graph_file_name = "./%s_stats_full.png" % (tail)
    print "out_graph_file_name: %s" % out_graph_file_name
    print "out_stats_graph_file_name: %s" % out_stats_graph_file_name
    out_data_file_name = "%s/%s/%s_full.json" % (head, tail, tail)
    print "out_data_file_name: %s" % out_data_file_name

    if os.path.exists(out_data_file_name):
        plot_data_file_name = out_data_file_name
        print "plot_data_file_name: %s" % plot_data_file_name

#    print "create connection\n"
#    db = MySQLdb.connect(
#        host="172.16.182.228",
#        user="tdesell",
#        passwd=passwd,
#        db="csg",
#        port=3306
#    )
#    print "connected\n"
#    cur = db.cursor()
#    cur.execute(
#        "\n".join((
#            "UPDATE gibbs_sampler",
#            "JOIN gibbs_post_status",
#            "   ON sampler_id = gibbs_sampler.id",
#            "SET run_ks = 'running',",
#            "run_ks_exception = NULL,",
#            "run_ks_job_id = %s",
#            "WHERE name=%s"
#        )),
#        (os.environ['PBS_JOBID'], tail)
#    )

    try:
        (file_values, total_files) = get_files_in_dir(in_dir)
        data = process_files_intrawalk(file_values, in_dir, out_data_file_name, total_files, plot_data_file_name)
        #plot_data(data, out_graph_file_name, title)
        plot_data_stats(data, out_stats_graph_file_name, title)
#        cur.execute(
#            "\n".join((
#                "UPDATE gibbs_sampler",
#                "JOIN gibbs_post_status",
#                "   ON sampler_id = gibbs_sampler.id",
#                "SET run_ks = 'rest',",
#                "run_ks_job_id = NULL",
#                "WHERE name=%s"
#            )),
#            (os.environ['PBS_JOBID'], tail)
#        )

    except Exception, ex:
#        cur.execute(
#            "\n".join((
#                "UPDATE gibbs_sampler",
#                "JOIN gibbs_post_status",
#                "   ON sampler_id = gibbs_sampler.id",
#                "SET run_ks = 'failed',",
#                "run_ks_job_id = NULL",
#                "run_ks_exception = %s",
#                "WHERE name = %s"
#            )),
#            (ex, tail)
#        )
        print "ERROR: %s" % ex

def ks_2(file_a=None, file_b=None, data_a=None, data_b=None):
    """Perform the Kolmogorov Smirnov 2 sample comparison"""
    start_time = time.time()
    if data_a is None:
        if file_a is None:
            print "file_a is None and data_a is None"
            raise ValueError()
        with open(file_a, 'r') as thing_a:
            #data_a = np.fromstring(thing_a.read(), dtype=int, sep=",").sort()

            data_a = np.fromstring(",".join([line.rstrip() for line in thing_a]), dtype=int, sep=",")
            print "data_a_shape: %s" % data_a.shape

    if data_b is None:
        if file_b is None:
            raise ValueError("file_b is None and data_b is None")
        with open(file_b, 'r') as thing_b:
            #data_b = np.fromstring(thing_b.read(), dtype=int, sep=",").sort()

            data_b = np.fromstring(",".join([line.rstrip() for line in thing_b]), dtype=int, sep=",")
            print "data_b_shape: %s" % data_b.shape

    data_read_time = time.time()
    dr_time = (data_read_time - start_time)
    ks_data = ks_2samp(data_a, data_b)
    ks_time = time.time() - data_read_time
    return [data_a, data_b, ks_data, dr_time, ks_time]

def get_files_in_dir(in_dir):
    """Organize files by walk and sort"""

    parents = sorted(os.listdir(in_dir))
    file_values = {}
    for parent in parents:
        try:
            files = sorted(os.listdir("%s/%s" % (in_dir, parent)))
        except Exception as ex:
            print ex

        pattern = re.compile(r"^samples_(\d+)_(\d+)$")

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

    for walk_name in file_values:
        file_values[walk_name] = sorted(file_values[walk_name])

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

    #(out_path, name) = os.path.split(in_dir)
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
                file_name_b = "_".join(["%s/walk_%s/samples_%s_%s" % (in_dir, walk_str, walk_str, step_str)])
                if walk_str not in data:
                    data[walk_str] = {}
                continue

            #print "data4: %s" % data

            file_name_a = file_name_b
            file_name_b = "_".join(["%s/walk_%s/samples_%s_%s" % (in_dir, walk_str, walk_str, step_str)])
            #print "step: %s" % step
            #print "str step: %s" % str(step)
            #print "step type: %s" % type(step)
            #print "data5: %s" % data
            #print "data[walk]: %s" % data[walk_str]
            #print "keys: %s" % data[walk_str].keys()
            #print "type keys: %s" % type(data[walk_str].keys())

            #XXX diasbled for test of 1st line bug
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

def plot_data(data, out_graph_file_name, title):
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

    plt.subplot(211)
    plt.title("%s Raw" % title)
    plt.savefig(out_graph_file_name)
    #plt.show()
    #print "KS results: %s" % data

def plot_data_stats(data, out_graph_file_name, title):
    """graph data reduction and error analysis"""


    max_step = 0
    max_walk = None
    for walk in data:
        if len(data[walk]) > max_step:
            max_step = len(data[walk])
            max_walk = walk
    if max_walk is None:
        print "max_walk is None, data: %s" % data
    x_data = sorted([int(n) for n in data[max_walk]])
    outer_y1_avg = np.empty([len(x_data)], dtype=float)
    outer_y1_std = np.empty([len(x_data)], dtype=float)
    outer_y1_min = np.empty([len(x_data)], dtype=float)
    outer_y1_max = np.empty([len(x_data)], dtype=float)

    outer_y2_avg = np.empty([len(x_data)], dtype=float)
    outer_y2_std = np.empty([len(x_data)], dtype=float)
    outer_y2_min = np.empty([len(x_data)], dtype=float)
    outer_y2_max = np.empty([len(x_data)], dtype=float)

    step_count = 0
    for step in x_data:
        str_step = str(step)

        y1_total = 0
        y2_total = 0
        y1_sos = 0
        y2_sos = 0
        count = 0
        y1_max = 0
        y2_max = 0
        y1_min = 100000
        y2_min = 100000


        for walk in data:
            if str_step in data[walk]:
                val_1 = data[walk][str_step][0][0]
                val_2 = data[walk][str_step][0][1]
                y1_total += val_1
                y2_total += val_2

                y1_sos += val_1 * val_1
                y2_sos += val_2 * val_2

                if y1_max < val_1:
                    y1_max = val_1
                if y2_max < val_2:
                    y2_max = val_2

                if val_1 < y1_min:
                    y1_min = val_1
                if val_2 < y2_min:
                    y2_min = val_2

                count += 1

        outer_y1_avg[step_count] = (y1_total / count)
        outer_y1_std[step_count] = (math.sqrt((y1_sos - y1_total * y1_total / count) / count))
        outer_y1_min[step_count] = (y1_min)
        outer_y1_max[step_count] = (y1_max)

        outer_y2_avg[step_count] = (y2_total / count)
        try:
            outer_y2_std[step_count] = (math.sqrt((y2_sos - y2_total * y2_total / count) / (count - 1)))
        except Exception, ex:
            print "%s" % ex
            print "sos - suare of sums: %s" % (y2_sos - y2_total * y2_total / count)
            print "(sqrt((%s - %s * %s / %s) / %s))" % (y2_sos, y2_total, y2_total, count, count)
            outer_y2_std[step_count] = 0
        outer_y2_min[step_count] = (y2_min)
        outer_y2_max[step_count] = (y2_max)
        step_count += 1

    fig_2 = plt.figure(2)
    #pivot
    #fig_2, (axes_1, axes_2) = plt.subplots(nrows=2, sharex=True)
    axes_1 = plt.subplot(211)
    patch = plt.Rectangle((0, 0), 0, 0, label="1 std dev range", color="black", alpha=0.5)
    axes_1.add_patch(patch)
    axes_1.set_ylim(ymin=0.0)
    plt.grid(True)
    plt.ylabel("Distance")
    #plt.plot(x_data, outer_y1_avg, "b", x_data, outer_y1_min, "ro", x_data, outer_y1_max, "go", x_data, outer_y1_avg - outer_y1_std, "r", x_data, outer_y1_avg + outer_y1_std, "g")
    plt.plot(x_data, outer_y1_avg, lw=3, label="Mean", color="black", ls="-")
    plt.fill_between(x_data, outer_y1_avg - outer_y1_std, outer_y1_avg + outer_y1_std, facecolor="black", alpha=0.25, label="1 deviation range")
    plt.plot(x_data, outer_y1_min, lw=2, label="Min", color="black", ls="-.")
    plt.plot(x_data, outer_y1_max, lw=2, label="Max", color="black", ls="--")
    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0., title=title)
    #plt.legend(bbox_to_anchor=[0.5, 0.5], loc='center', ncol=4, mode="expand", borderaxespad=0.)
    #fig_2.legend(bbox_to_anchor=[0.5, 0.5], loc='center', ncol=4)
    axes_1.legend(bbox_to_anchor=[0.5, 0], loc='upper center', ncol=4, borderaxespad=0.25, fancybox=True, fontsize='small')
    axes_1.set_xticklabels([])
    axes_1.set_yscale("log")
    axes_1.set_adjustable("datalim")
    axes_1.set_ylim(1e-3, 1e0)
    #axes_1.set_aspect(1)
    #plt.semilogy(x_data, outer_y1_avg, "b", x_data, outer_y1_min, "ro", x_data, outer_y1_max, "go", x_data, outer_y1_avg - outer_y1_std, "r", x_data, outer_y1_avg + outer_y1_std, "g")
    axes_2 = plt.subplot(212)
    #axes_2.set_yscale("log")
    #axes_2.set_adjustable("datalim")
    #axes_2.set_ylim(1e-2, 1e0)
    #axes_2.set_aspect(1)

    axes_2.set_ylim(ymin=0.0, ymax=1.1)
    plt.grid(True)
    plt.xlabel("Step")
    plt.ylabel("Probability")
    #plt.plot(x_data, outer_y2_avg, "b", x_data, outer_y2_min, "ro", x_data, outer_y2_max, "go", x_data, outer_y2_avg - outer_y2_std, "r", x_data, outer_y2_avg + outer_y2_std, "g")
    #plt.semilogy(x_data, outer_y2_avg, "b", x_data, outer_y2_min, "ro", x_data, outer_y2_max, "go", x_data, outer_y2_avg - outer_y2_std, "r", x_data, outer_y2_avg + outer_y2_std, "g")

    plt.plot(x_data, outer_y2_avg, lw=3, label="Mean", color="black", ls="-")
    plt.fill_between(x_data, outer_y2_avg - outer_y2_std, outer_y2_avg + outer_y2_std, facecolor="black", alpha=0.25, label="1 deviation range")
    plt.plot(x_data, outer_y2_min, lw=2, label="Min", color="black", ls="-.")
    plt.plot(x_data, outer_y2_max, lw=2, label="Max", color="black", ls="--")
    plt.subplot(211)
    #plt.title(os.path.splitext(os.path.basename(out_graph_file_name))[0])
    plt.title(" ".join(title.split("_")))
    plt.savefig("%s" % out_graph_file_name)



if __name__ == "__main__":
    main()

