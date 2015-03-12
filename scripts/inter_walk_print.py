#!/usr/bin/python
""" sum and possible average the motifs we've found.  consider weighted average"""

import json
import re
from optparse import OptionParser
import os


#no match?     0, 0     988 atgtt CCTCCT tttcc     993 0.2110 > 'PEF1 chr1 32109718 32110718'
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


    for walk in files:
    #for walk in result:

        with open("%s/%s" % (options.sample_dir, walk), "r") as in_file:
            try:
                gen = (process_line(motif_stats, line) for line in in_file)

                while 1:
                    next(gen)

            except StopIteration:
                print "file %s done" % walk
                in_file.close()


    completion(motif_stats)
    with open("%s/motif_stats_%s_%s.json" % (options.sample_dir, tail, options.step), "w") as out_file:
        json.dump(motif_stats, out_file, indent=4, separators=(',', ': '))
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
#XXX we can process for neigbors or other stats later
def process_line(motif_stats, line):
    """ read a line, process it and add it to the motifs"""

    if HEADER_R.match(line):
        return
    else:
        body_m = BODY_R.match(line)
        if body_m is not None:
            (seq_num, motif_num, start, lneighbor, motif, rneighbor, end, pct, gene_id) = body_m.groups()

            seq_num = int(seq_num)
            motif_num = int(motif_num)

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

            if not (leboxc or leboxg):
                return

            pct = float(pct) * 100
            if gene_id not in motif_stats["by_gene"]:
                motif_stats["by_gene"][gene_id] = {}

            if motif not in motif_stats["by_gene"][gene_id]:
                motif_stats["by_gene"][gene_id][super_motif] = {
                    "sequence" : seq_num,
                    "total_pct" : 0.0,
                    "count" : 0.0,
                    "max_pct" : 0.0,
                    "min_pct" : 100.0,
                    "sos_pct" : 0.0,
                    "max_start" : None,
                    "eboxc" : leboxc,
                    "eboxg" : leboxg
                }

            lmod = motif_stats["by_gene"][gene_id][super_motif]
            lmod["total_pct"] = lmod["total_pct"] + pct
            lmod["count"] = lmod["count"] + 1

            if pct > lmod["max_pct"]:
                lmod["max_pct"] = pct
                lmod["max_start"] = start

            if pct < lmod["min_pct"]:
                lmod["min_pct"] = pct

            lmod["sos_pct"] = lmod["sos_pct"] + pct * pct

            if super_motif not in motif_stats["global"]:
                motif_stats["global"][super_motif] = {
                    "total_pct" : 0.0,
                    "count" : 0,
                    "max_pct" : 0.0,
                    "min_pct" : 100.0,
                    "sos_pct" : 0.0,
                    "max_gene_id" : None,
                    "max_gene_start" : None,
                    "eboxc" : leboxc,
                    "eboxg" : leboxg
                }

            gmod = motif_stats["global"][super_motif]
            gmod["total_pct"] = gmod["total_pct"] + pct
            gmod["count"] = gmod["count"] + 1

            if pct > gmod["max_pct"]:
                gmod["max_pct"] = pct
                gmod["max_gene_id"] = gene_id
                gmod["max_start"] = start

            if pct < gmod["min_pct"]:
                gmod["min_pct"] = pct

            gmod["sos_pct"] = gmod["sos_pct"] + pct * pct

        else:
            print "no match? %s" % line

def completion(motif_stats):
    """ finalize statistics"""
    for gene_id in motif_stats["by_gene"]:
        for super_motif in motif_stats["by_gene"][gene_id]:

            lmotif = motif_stats["by_gene"][gene_id][super_motif]
            lmotif["avg_pct_best"] = lmotif["total_pct"] / lmotif["count"]
            lmotif["avg_pct"] = lmotif["total_pct"] / motif_stats["sample_count"]
            lmotif["std_dev_best"] = (lmotif["sos_pct"] - lmotif["total_pct"] * lmotif["total_pct"]) / lmotif["count"]
            lmotif["std_dev"] = (
                (lmotif["sos_pct"] - lmotif["total_pct"] * lmotif["total_pct"])
                / (motif_stats["sample_count"] * motif_stats["motif_num"])
            )

            lmotif["eboxc"] = False
            lmotif["eboxg"] = False
            if EBOXC_R.match(super_motif):
                lmotif["eboxc"] = True

            if EBOXG_R.match(super_motif):
                lmotif["eboxg"] = True

    for super_motif in motif_stats["global"]:

        gmotif = motif_stats["global"][super_motif]
        gmotif["avg_pct_best"] = gmotif["total_pct"] / gmotif["count"]
        gmotif["avg_pct"] = gmotif["total_pct"] / motif_stats["sample_count"]
        gmotif["std_dev_best"] = (gmotif["sos_pct"] - gmotif["total_pct"] * gmotif["total_pct"]) / gmotif["count"]
        gmotif["std_dev"] = (
            (gmotif["sos_pct"] - gmotif["total_pct"] * gmotif["total_pct"])
            / (motif_stats["sample_count"] * motif_stats["motif_num"])
        )
        gmotif["eboxc"] = False
        gmotif["eboxg"] = False
        if EBOXC_R.match(super_motif):
            gmotif["eboxc"] = True

        if EBOXG_R.match(super_motif):
            gmotif["eboxg"] = True


if __name__ == "__main__":
    main()

