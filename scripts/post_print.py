#!/usr/bin/python
""" sum and possible average the motifs we've found.  consider weighted average"""

import json
import csv
import re
from optparse import OptionParser
import os
import math

#no match?     0, 0     988 atgtt CCTCCT tttcc     993 0.2110 > 'ISG15 chr1 948346 949346 948846 +'
HEADER_R = re.compile("^[FR]|$")
BODY_R = re.compile(r"^\s+(\d+),\s*(\d+)\s+(\d+)\s+([a-z]*)\s*([A-Z]+)\s*([a-z]*)\s+(\d+)\s+([0-9\.]+)\s+>\s+'([^']+)'$")

#these guys are the proof in the pudding
EBOXC_R = re.compile(r"CACCTG", flags=re.IGNORECASE)
EBOXG_R = re.compile(r"CAGGTG", flags=re.IGNORECASE)

def main():
    """ cleanup"""
    #usage = "usage: %prog [options] arg"
    usage = ""
    parser = OptionParser(usage)
    parser.add_option("--sample_dir", dest="sample_dir", help="What sample directory should be used?", metavar="STRING", type="string")
    parser.add_option("--best_pct", dest="best_pct", help="What was the pct cutoff used to generate the data?", metavar="FLOAT", type="float")
    parser.add_option("--ebox_only", dest="ebox_only", help="Only return motifs containing ebox with motif or neighbors", action="store_true", default=False)

    (options, args) = parser.parse_args()
    if not options.sample_dir or not options.best_pct:
        parser.print_help()
        parser.error("missing options: %s" % options)

    #motif_files = ["./walk_200010_steps_30000.motifs"]
    #for motif_file in motif_files:
    search = None
    search = re.compile("^walk_steps_all_(\d+_\d+).motifs")

    files = [f for f in os.listdir(options.sample_dir) if re.match(search, f)]

    (head, tail) = os.path.split(options.sample_dir)
    if tail == '':
        (head, tail) = os.path.split(head)



    # 23, 1      31 tgtgcTCTCCTaaacc      36 0.1702 > 'PNP chr14 20937037 20938037'
    motif_stats = {
        "by_gene" : {},
        "best_pct" : options.best_pct,
        "source" : tail,
    }
    print files

    #for walk in files:

    print "input: %s" % files[0]
    with open("%s/%s" % (options.sample_dir, files[0]), "r") as in_file:
        try:
            gen = (process_line(motif_stats, line, options) for line in in_file)

            while 1:
                next(gen)

        except StopIteration:
            print "file %s done" % files[0]
            in_file.close()

    #print "motif_stats: %s" % motif_stats
    out = completion(motif_stats, options)
    match = re.match(search, files[0])

    out_name = "./motif_stats_%s_%s.csv" % (tail, match.group(1))
    print "output written to: %s" % out_name

    with open(out_name, "w") as out_file:
        #json.dump(motif_stats, out_file, indent=4, separators=(',', ': '))
        writer = csv.writer(out_file, dialect=csv.excel_tab)
        for row in out:
            writer.writerow(row.split("\t"))
        out_file.close()

#motif_stats plan
#{
#  "best_pct" : .01,
#  "by_gene" : {
#    $gene_id : {
#      $location : {
#        total_pct : 1234,
#        count : 500,
#        max_pct : 12,
#        min_pct : 12,
#        sos_pct : 1234,
#        max_start : 11111
#      },
#    },
#  },
#}
# we can process for neigbors or other stats later
def process_line(motif_stats, line, options):
    """ read a line, process it and add it to the motifs"""

    if HEADER_R.match(line):
        return
    else:
        body_m = BODY_R.match(line)
        if body_m is not None:
            (seq_num, motif_num, start, lneighbor, motif, rneighbor, end, pct, gene_id) = body_m.groups()

            seq_num = int(seq_num)

            #combine neighbors with motif for sorting and searching for ebox
            super_motif = "%s%s%s" % (lneighbor, motif, rneighbor)

            leboxc = False
            leboxg = False
            if EBOXC_R.match(super_motif):
                leboxc = True

            if EBOXG_R.match(super_motif):
                leboxg = True

            if options.ebox_only and not (leboxc or leboxg):
                return

            pct = float(pct)

            # so the gene_id is everything in the '' here
            #no match?     0, 0     988 atgtt CCTCCT tttcc     993 0.2110 > 'ISG15 chr1 948346 949346 948846 +'
            #this is not so efficient, we only really need ISG15.  It won't appear again right?
            #we do need all of this information but a shorter hash key would be good
            #what is in the gene_id:  gene symbol, chromosome, region start, region end, transcription start, strand
            #note the gene_symbol and transcription start and strand may be a comma spearated list
            #if you have to work with gene symbols that already contain punctuation, the guy that chose the symbol
            #needs a talking to

            (gene_symbol, chrom, region_start, region_end, tss, strand) = gene_id.split(" ")
            if gene_id not in motif_stats["by_gene"]:
                motif_stats["by_gene"][gene_id] = {}
            location = str(int(start) + int(region_start))

            if location not in motif_stats["by_gene"][gene_id]:

                tss_offsets = []
                for single_tss in tss.split(","):
                    #start is an offset from region_start
                    #to change it to an offset from tss
                    #we need to find the difference of tss and region_start
                    tss_offsets.append(str(int(start) - (int(single_tss) - int(region_start))))

                motif_stats["by_gene"][gene_id][location] = {
                    "total_pct" : 0.0,
                    "count" : 0.0,
                    "max_pct" : pct,
                    "min_pct" : pct,
                    "cacctg" : leboxc,
                    "caggtg" : leboxg,
                    "start" : start,
                    "location" : location,
                    "sos_pct" : 0.0,
                    "lneighbor" : lneighbor,
                    "rneighbor" : rneighbor,
                    "motif" : motif.lower(),
                    "gene_symbol" : gene_symbol,
                    "chrom" : chrom,
                    "region_start" : region_start,
                    "region_end" : region_end,
                    "tss" : tss,
                    "strand" : strand,
                    "tss_offset" : ",".join(tss_offsets)
                }

            lmod = motif_stats["by_gene"][gene_id][location]
            lmod["total_pct"] = lmod["total_pct"] + pct
            lmod["count"] = lmod["count"] + 1

            if pct > lmod["max_pct"]:
                lmod["max_pct"] = pct

            if pct < lmod["min_pct"]:
                lmod["min_pct"] = pct

            lmod["sos_pct"] = lmod["sos_pct"] + pct * pct

        else:
            print "no match? %s" % line

def completion(motif_stats, options):
    """ finalize statistics"""
    print "motif_stats: %s" % motif_stats
    readable = []
    for gene_id in motif_stats["by_gene"]:
        for location in motif_stats["by_gene"][gene_id]:

            lmotif = motif_stats["by_gene"][gene_id][location]
            lmotif["avg_pct_best"] = lmotif["total_pct"] / lmotif["count"]
#            lmotif["avg_pct"] = lmotif["total_pct"] / motif_stats["sample_count"]
            lmotif["std_dev_best"] = math.sqrt((lmotif["sos_pct"] - lmotif["total_pct"] * lmotif["total_pct"] / lmotif["count"])/lmotif["count"])
            #if lmotif["count"] > 100 and lmotif["avg_pct_best"] > 10:
            if lmotif["avg_pct_best"] > options.best_pct:
                readable.append("\t".join((
                    lmotif["chrom"],
                    lmotif["gene_symbol"],
                    lmotif["tss"],
                    lmotif["motif"],
                    lmotif["location"],
                    "%.2f" % (lmotif["avg_pct_best"]),
                    lmotif["region_start"],
                    lmotif["region_end"],
                    lmotif["tss_offset"],
                    lmotif["lneighbor"],
                    lmotif["rneighbor"],
                    "%s%s%s" % (lmotif["lneighbor"], lmotif["motif"], lmotif["rneighbor"]),
                    str(lmotif["cacctg"]),
                    str(lmotif["caggtg"]),
                    lmotif["strand"],
                    "%d" % (lmotif["count"])
                )))
    print "presort: %s" % readable
    sort_nicely(readable)
    print "postsort: %s" % readable
    readable.insert(0, ("\t".join((
        "chromosome",
        "gene_symbol",
        "tss",
        "motif",
        "location",
        "avg_pct",
        "interval_start",
        "interval_end",
        "tss_offset",
        "left_neighbor",
        "right_neighbor",
        "combined_motif",
        "cacctg",
        "caggtg",
        "strand",
        "count"
    ))))
    return readable

def tryint(string):
    """is it a int?"""
    try:
        return int(string)
    except:
        return string

def alphanum_key(string):
    """find ints strings"""
    return [tryint(char) for char in re.split('([0-9]+)', string)]

def sort_nicely(unsorted_list):
    """natural sort"""
    unsorted_list.sort(key=alphanum_key)

if __name__ == "__main__":
    main()
