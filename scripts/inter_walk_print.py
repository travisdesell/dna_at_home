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
    parser.add_option("--step", dest="step", help="Which step should be compared?", metavar="INT", type="int")
    parser.add_option("--sample_dir", dest="sample_dir", help="What sample directory should be used?", metavar="STRING", type="string")
    parser.add_option("--best_pct", dest="best_pct", help="What was the pct cutoff used to generate the data?", metavar="FLOAT", type="float")
    parser.add_option("--ebox_only", dest="ebox_only", help="Only return motifs containing ebox with motif or neighbors", action="store_true", default=False)
    parser.add_option("--tabbed_csv", dest="tabbed_csv", help="Produce a tab delimited csv instead of latex format output", action="store_true", default=False)

    (options, args) = parser.parse_args()
    if not options.step or not options.sample_dir or not options.best_pct:
        parser.print_help()
        parser.error("missing options: %s" % options)

    #motif_files = ["./walk_200010_steps_30000.motifs"]
    #for motif_file in motif_files:
    search = re.compile("^.*_%s.motifs" % options.step)
    files = [f for f in os.listdir(options.sample_dir) if re.match(search, f)]

    (head, tail) = os.path.split(options.sample_dir)
    if tail == '':
        (head, tail) = os.path.split(head)



    # 23, 1      31 tgtgcTCTCCTaaacc      36 0.1702 > 'PNP chr14 20937037 20938037'
    motif_stats = {
        "by_gene" : {},
        "global" : {},
        "sample_count" : len(files),
        "best_pct" : options.best_pct,
        "source" : tail,
        "step" : options.step,
        "motif_num" : -1,
        "seq_num" : -1
    }
    print files

    for walk in files:

        print "walk: %s" % walk
        with open("%s/%s" % (options.sample_dir, walk), "r") as in_file:
            try:
                gen = (process_line(motif_stats, line, options) for line in in_file)

                while 1:
                    next(gen)

            except StopIteration:
                print "file %s done" % walk
                in_file.close()

    #print "motif_stats: %s" % motif_stats
    out = completion(motif_stats, options)
    #print "out: %s" % out
    out_name = "%s/motif_stats_%s_%s.json"
    if options.tabbed_csv:
        out_name = "%s/motif_stats_%s_%s.csv"

    print "output written to: %s" % out_name
    with open(out_name % (options.sample_dir, tail, options.step), "w") as out_file:
        #json.dump(motif_stats, out_file, indent=4, separators=(',', ': '))
        if options.tabbed_csv:
            writer = csv.writer(out_file, dialect=csv.excel_tab)
            for row in out:
                writer.writerow(row.split("\t"))
        else:
            json.dump(out, out_file, indent=4, separators=(',', ': '))
        out_file.close()

#motif_stats plan
#{
#  "sample_count" : 1000,
#  "best_pct" : .01,
#  "by_gene" : {
#    $gene_id : {
#      $motif : {
#        total_pct : 1234,
#        count : 500,
#        max_pct : 12,
#        min_pct : 12,
#        sos_pct : 1234,
#        max_start : 11111
#      },
#    },
#  },
#  "global" : {
#    $motif : {
#      total_pct : 1234,
#      count : 500,
#      max_pct : 12,
#      min_pct : 12,
#      sos_pct : 1234,
#      max_gene_id : asdf,
#      max_gene_start
#    },
#  }
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
            motif_num = int(motif_num)

            #combine neighbors with motif for sorting and searching for ebox
            super_motif = "%s%s%s" % (lneighbor, motif, rneighbor)


            if motif_num > motif_stats["motif_num"]:
                motif_stats["motif_num"] = motif_num + 1

            if seq_num > motif_stats["seq_num"]:
                motif_stats["seq_num"] = seq_num + 1

            leboxc = False
            leboxg = False
            if EBOXC_R.match(super_motif):
                leboxc = True

            if EBOXG_R.match(super_motif):
                leboxg = True

            if options.ebox_only and not (leboxc or leboxg):
                return

            pct = float(pct) * 100

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

            if super_motif not in motif_stats["by_gene"][gene_id]:

                tss_offsets = []
                for single_tss in tss.split(","):
                    
                    tss_offsets.append(str(int(single_tss) - (int(region_start) + int(start))))

                motif_stats["by_gene"][gene_id][super_motif] = {
                    "total_pct" : 0.0,
                    "count" : 0.0,
                    "max_pct" : 0.0,
                    "min_pct" : 100.0,
                    "cacctg" : leboxc,
                    "caggtg" : leboxg,
                    "start" : start,
                    "sos_pct" : 0.0,
                    "lneighbor" : lneighbor,
                    "rneighbor" : rneighbor,
                    "motif" : motif,
                    "gene_symbol" : gene_symbol,
                    "chrom" : chrom,
                    "region_start" : region_start,
                    "region_end" : region_end,
                    "tss" : tss,
                    "strand" : strand,
                    "tss_offset" : ",".join(tss_offsets)
                }
#
            lmod = motif_stats["by_gene"][gene_id][super_motif]
            lmod["total_pct"] = lmod["total_pct"] + pct
            lmod["count"] = lmod["count"] + 1
#
            if pct > lmod["max_pct"]:
                lmod["max_pct"] = pct
#                lmod["max_start"] = start
#
            if pct < lmod["min_pct"]:
                lmod["min_pct"] = pct
#
            lmod["sos_pct"] = lmod["sos_pct"] + pct * pct
#
#            if super_motif not in motif_stats["global"]:
#                motif_stats["global"][super_motif] = {
#                    "total_pct" : 0.0,
#                    "count" : 0,
#                    "max_pct" : 0.0,
#                    "min_pct" : 100.0,
##                    "sos_pct" : 0.0,
#                    "eboxc" : leboxc,
#                    "eboxg" : leboxg,
#                    "genes" : {}
#                }
#
#
#            gmod = motif_stats["global"][super_motif]
#            gmod["total_pct"] = gmod["total_pct"] + pct
#            gmod["count"] = gmod["count"] + 1
#
#            if pct > gmod["max_pct"]:
#                gmod["max_pct"] = pct
# #               gmod["max_gene_id"] = gene_id
# #               gmod["max_start"] = start
#
#            if pct < gmod["min_pct"]:
#                gmod["min_pct"] = pct
#
#            if gene_id not in gmod["genes"]:
#                gmod["genes"][gene_id] = line
##            gmod["sos_pct"] = gmod["sos_pct"] + pct * pct
#
        else:
            print "no match? %s" % line

def completion(motif_stats, options):
    """ finalize statistics"""
    readable = []
    for gene_id in motif_stats["by_gene"]:
        for super_motif in motif_stats["by_gene"][gene_id]:

            lmotif = motif_stats["by_gene"][gene_id][super_motif]
            lmotif["avg_pct_best"] = lmotif["total_pct"] / lmotif["count"]
#            lmotif["avg_pct"] = lmotif["total_pct"] / motif_stats["sample_count"]
            lmotif["std_dev_best"] = math.sqrt((lmotif["sos_pct"] - lmotif["total_pct"] * lmotif["total_pct"] / lmotif["count"])/lmotif["count"])
            if lmotif["count"] > 100 and lmotif["avg_pct_best"] > 10:
                if options.tabbed_csv:
                    #estimated tss offset is due to combined genes.  if we had 2 tss...
                    readable.append("\t".join((
                        "%d" % (lmotif["count"]),
                        "%.2f" % (lmotif["avg_pct_best"]),
                        lmotif["gene_symbol"],
                        lmotif["chrom"],
                        lmotif["region_start"],
                        lmotif["region_end"],
                        lmotif["tss"],
                        lmotif["tss_offset"],
                        lmotif["lneighbor"],
                        lmotif["motif"],
                        lmotif["rneighbor"],
                        super_motif,
                        str(lmotif["cacctg"]),
                        str(lmotif["caggtg"]),
                        lmotif["strand"]
                    )))
                else:
                    readable.append("& ".join((
                        "%d" % (lmotif["count"]),
                        "%.2f" % (lmotif["avg_pct_best"]),
                        gene_id,
                        str(lmotif["start"]),
                        super_motif,
                        str(lmotif["cacctg"]),
                        str(lmotif["caggtg"])
                    )))
    sort_nicely(readable)
    if options.tabbed_csv:
        readable.insert(0, ("\t".join((
            "count",
            "avg_pct",
            "gene_symbol",
            "chromosome",
            "interval_start",
            "interval_end",
            "tss",
            "tss_offset",
            "left_neighbor",
            "motif",
            "right_neighbor",
            "combined_motif",
            "cacctg",
            "caggtg",
            "strand"
        ))))
    return readable
#            lmotif["std_dev"] = ((
#                lmotif["sos_pct"] - lmotif["total_pct"] * lmotif["total_pct"]
#                / (motif_stats["sample_count"] * motif_stats["motif_num"])
#            )/(motif_stats["sample_count"] * motif_stats["motif_num"]))

#            lmotif["eboxc"] = False
#            lmotif["eboxg"] = False
#            if EBOXC_R.match(super_motif):
#                lmotif["eboxc"] = True

#            if EBOXG_R.match(super_motif):
#                lmotif["eboxg"] = True

#    for super_motif in motif_stats["global"]:
#
#        gmotif = motif_stats["global"][super_motif]
#        gmotif["avg_pct_best"] = gmotif["total_pct"] / gmotif["count"]
#        gmotif["avg_pct"] = gmotif["total_pct"] / motif_stats["sample_count"]
#        gmotif["std_dev_best"] = (gmotif["sos_pct"] - gmotif["total_pct"] * gmotif["total_pct"]) / gmotif["count"]
#        gmotif["std_dev"] = (
#            (gmotif["sos_pct"] - gmotif["total_pct"] * gmotif["total_pct"])
#            / (motif_stats["sample_count"] * motif_stats["motif_num"])
#        )
#        gmotif["eboxc"] = False
#        gmotif["eboxg"] = False
#        if EBOXC_R.match(super_motif):
#            gmotif["eboxc"] = True
#
#        if EBOXG_R.match(super_motif):
#            gmotif["eboxg"] = True
#
def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    return [tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    l.sort(key=alphanum_key)

if __name__ == "__main__":
    main()
