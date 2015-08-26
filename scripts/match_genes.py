#!/usr/bin/python
"""Match gene info to fastas"""

import re
from collections import defaultdict
from optparse import OptionParser


#every line from the cpp_result input file is a new gene site
class GeneSite:
    """a gene region datastructure"""
    def __init__(self, line, db_info):
        #cpp_result is tab separated
        #it contians the following fields: db.desc1, db.desc2, db.chr, db.start, db.end, q.desc1, q.hits, q.chr, q.start, q.end, match
        #that's
        #db.desc1, a gene id
        #db.desc2, a gene symbol:transcription start site
        #db.chr, the chromosome
        #dv.start, starting point of the gene region
        #db.end, ending point of the gene region
        #q.desc1, the cpp query, arbitrary identifier
        #q.hits, an expression rating from cppmatch
        #q.chr, the query chromsomoe same as before
        #q.start, a bin starting point, unused here
        #q.end, a bin ending point, unused here

        #since we really only care about the db portion of this, if we have a db_file
        #that's been filtered to the set of desired genes and we don't care about limiting
        #then we can just pass that in
        #split without addressing words > 4 is already effectively truncating the cpp results
        words = re.split("[ \t]", line)
        #separate the gene symbol and transcription start site
        (symbol, tss) = words[1].split(":")
        self.gene = symbol
        self.chromosome = words[2]
        self.start = int(words[3])
        self.end = int(words[4])

        self.tss = tss
        if db_info is not None:
            self.strand = db_info[symbol]
        elif words[5] == "+" or words[5] == "-":
            self.strand = words[5]
        else:
            print "we lost strand info"

        #determine a strand based on the direction of the region
        if self.start > self.end:
            #adjust for the direction of the region, gibbs sampler wants everything in the same direction
            #it will guess at directionality
            print "need to swap start/end for: %s" % self.gene
            tmp = self.start
            self.start = self.end
            self.end = tmp

    def __repr__(self):
        return "<%s %s %s %s %s %s>" % (self.gene, self.chromosome, self.start, self.end, self.tss, self.strand)

    def __str__(self):
        return "'%s %s %s %s %s %s'" % (self.gene, self.chromosome, self.start, self.end, self.tss, self.strand)

def main():
    """ let's get some options and determine how to run this thing"""

    usage = ""
    parser = OptionParser(usage)
    parser.add_option("--sites_file", dest="sites_file", help="What db file or cpp_result file to use?", metavar="STRING", type="string")
    parser.add_option("--chrom_dir", dest="chrom_dir", help="What directory has the chromosome fastas?", metavar="STRING", type="string")
    parser.add_option("--out", dest="out", help="What is the name of the output file?", metavar="STRING", type="string")
    parser.add_option("--limit", dest="limit", help="Should a limit be used?", metavar="INT", type="int")
    parser.add_option("--db_file", dest="db_file", help="Optional for recovering strand information if cpp file was used?", metavar="STRING", type="string")

    (options, args) = parser.parse_args()
    if not options.sites_file or not options.chrom_dir or not options.out:
        parser.print_help()
        parser.error("options: %s missing one of (sites_file, chrom_dir, out)" % options)

    db_info = None
    if options.db_file:
        db_info = get_db(options)
    sites = get_sites(options, db_info)
    (unduplicated_sites, chromosomes) = unduplicate(sites)
    chromosome_files = read_chrom(chromosomes, options)
    write_out(unduplicated_sites, chromosome_files, options)

def get_db(options):
    """create a kvp for symbol : strand"""
    db_info = {}
    with open(options.db_file) as db_file:
        for line in db_file:
            if line[0] == '#':
                continue

            line = line.strip()
            (id, symbol_tss, chrom, start, end, strand) = line.split("\t")
            (symbol, tss) = symbol_tss.split(":")
            db_info[symbol] = strand
    return db_info

def get_sites(options, db_info):
    """ create sites"""
    sites = []
    #open a cpp_result file
    with open(options.sites_file) as sites_file:
        #we may have arbitrarily limited results to the first X genes
        print "limit: %s" % options.limit
        next(sites_file)
        for line in sites_file:
            print "limit: %s, sites: %s\n" % (options.limit, len(sites))
            if options.limit is not None and len(sites) > options.limit:
                break
            #skip the header of a cpp_result
            if line[0] == "db.desc1":
                continue
            #skip the header of a db_file
            if line[0] == '#':
                continue

            sites.append(GeneSite(line, db_info))

    #expression rating was used to limit results if limits were desired
    #now we want to sort for our own purposes as we decide how to deal with overlapping genes
    sites = sorted(sites, key=lambda site: site.start)
    sites = sorted(sites, key=lambda site: site.chromosome)

    for site in sites:
        print site

    return sites

def read_chrom(chromosomes, options):
    """Read chromosomes in"""
    print "reading chromosome files."

    chromosome_files = defaultdict(str)
    #fasta data is separated in chromosome files
    #if size complexity were a problem, this would be changed to only have one chromosome in memory at a time
    for chromosome in chromosomes:

        with open(options.chrom_dir + chromosome + ".fa", 'r') as chromosome_file:

            print "read chromosome: %s" % chromosome
            contents = chromosome_file.read()
            pos = contents.index('\n')
            contents = contents[pos+1:]
            contents = re.sub('[\n\r \t]', '', contents)
            chromosome_files[chromosome] = contents


    print "chromosome files read"
    return chromosome_files

def unduplicate(sites):
    """remove duplicate sites and determine chromosomes needed"""
    unduplicated_sites = []
    chromosomes = []
    overlaps = 0
    pos = 0
    while pos < len(sites):
    #    print sites[pos]
        if sites[pos].chromosome not in chromosomes:
            chromosomes.append(sites[pos].chromosome)

        pos_start = pos
        pos_end = pos+1
        #combine the overlapped sites
        while pos_end < len(sites) and sites[pos_start].chromosome == sites[pos_end].chromosome and sites[pos_start].end > sites[pos_end].start:
            overlaps += 1
            #concatenate to the gene field
            sites[pos_start].gene += "," + sites[pos_end].gene
            #concatenate to the tss field
            sites[pos_start].tss += "," + sites[pos_end].tss
            sites[pos_start].end = sites[pos_end].end
            sites[pos_start].strand += "," + sites[pos_end].strand
            pos_end += 1

        unduplicated_sites.append(sites[pos])
        pos = pos_end

    print len(sites), "total gene sites."
    print overlaps, "overlaps."
    print len(unduplicated_sites), "unduplicated gene sites."

    return (unduplicated_sites, chromosomes)

def write_out(unduplicated_sites, chromosome_files, options):
    """put it all together now"""
    #recover the fasta data
    with open(options.out, 'w') as out:
        for unduplicated_site in unduplicated_sites:
            #contents_len = len(chromosome_files[unduplicated_site.chromosome])
            #print "chromosome length: %d, gene start: %d, gene end: %d" % (contents_len, unduplicated_site.start, unduplicated_site.end)

            snip = chromosome_files[unduplicated_site.chromosome][unduplicated_site.start:unduplicated_site.end]

            if snip == '':
                print "ERROR: snip was '' for: ", unduplicated_site
            else:
                out.write("> %s\n" % unduplicated_site)
                out.write("%s\n" % snip)
                out.write("\n")


if __name__ == "__main__":
    main()
