#!/usr/bin/python
"""Match gene info to fastas"""

from optparse import OptionParser
#import matplotlib.pyplot as plt
#import numpy as np
import heatmap

#every line from the cpp_result input file is a new gene site
class GeneSite:
    """a gene region datastructure"""
    def __init__(self, line):
        (
            self.count,
            self.avg_pct,
            self.gene_symbol,
            self.chromosome,
            self.interval_start,
            self.interval_end,
            self.tss,
            self.tss_offset,
            self.left_neighbor,
            self.motif,
            self.right_neighbor,
            self.combined_motif,
            self.cacctg,
            self.caggtg,
            self.strand
        ) = line.split("\t")

    def __repr__(self):
        return "<%s %s %s %s %s %s>" % (self.gene_symbol, self.chromosome, self.interval_start, self.interval_end, self.tss, self.strand)

    def __str__(self):
        return "'%s %s %s %s %s %s'" % (self.gene_symbol, self.chromosome, self.interval_start, self.interval_end, self.tss, self.strand)



def main():
    """ let's get some options and determine how to run this thing"""

    usage = ""
    parser = OptionParser(usage)
    parser.add_option("--in_file", dest="in_file", help="What data file to use?", metavar="STRING", type="string")
#    parser.add_option("--chrom_dir", dest="chrom_dir", help="What directory has the chromosome fastas?", metavar="STRING", type="string")
    parser.add_option("--out", dest="out", help="What is the name of the output file?", metavar="STRING", type="string")
    parser.add_option("--dot_size", dest="dot_size", help="dot_size", metavar="INT", type="int", default=15)
    parser.add_option("--opacity", dest="opacity", help="opacity", metavar="INT", type="int", default=128)
#    parser.add_option("--db_file", dest="db_file", help="Optional for recovering strand information if cpp file was used?", metavar="STRING", type="string")

    (options, args) = parser.parse_args()
    if not options.in_file or not options.out:
        parser.print_help()
        parser.error("options: %s missing one of (in_file, out)" % options)

    sites = parse_file(options.in_file)

    coords = gen_coordinates(sites)
    print coords


    hm = heatmap.Heatmap()
    img = hm.heatmap(coords, dotsize=options.dot_size, opacity=options.opacity)
    img.save("%s.png" % options.out)

#    x = []
#    y = []
#    for coord in coords:
#        x.append(coord[0])
#        y.append(coord[1])
#
#    heatmap, xedges, yedges = np.histogram2d(x, y, bins=50)
#    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#
#
#    plt.clf()
#    plt.imshow(heatmap, extent=extent)
#    plt.show()

def gen_coordinates(sites):
    """generate a list of tuples of coordinates for the heatmap"""
    coords = []

    for count in range(len(sites)):
        site = sites[count]
        print "site: %s" % site
        print "count: %s, tss_offset: %s" % (count, site.tss_offset)
        #coords.append((int(count), int(site.tss_offset)))
        coords.append((int(site.tss_offset), count))

    return coords

def parse_file(in_file):
    """split the file"""

    sites = []
    with open(in_file, "r") as in_data:
        next(in_data)
        for line in in_data:
            for site in parse_line(line.rstrip()):
                sites.append(site)

    return sites

def parse_line(line):
    """split the line or handle multigene lines"""

    sites = []
    if "," in line:
        print "more than one gene on this line"
        lines = split_line(line)
        for subline in lines:
            sites.append(GeneSite(subline))
    else:
        print "do something to a line"

        sites.append(GeneSite(line))

    return sites


def split_line(line):
    """handle multigene"""
    (
        count,
        avg_pct,
        gene_symbol,
        chromosome,
        interval_start,
        interval_end,
        tss,
        tss_offset,
        left_neighbor,
        motif,
        right_neighbor,
        combined_motif,
        cacctg,
        caggtg,
        strand
    ) = line.split("\t")

    gene_symbols = gene_symbol.split(",")
    tsses = tss.split(",")
    tss_offsets = tss_offset.split(",")
    strands = strand.split(",")

    lines = []
    for count in range(len(gene_symbols)):
        print "count: %s" % count
        print "gene_symbols: %s" % gene_symbols
        print "gene_symbols[count]: %s" % gene_symbols[count]
        print "strands: %s" % strands
        print "strands[count]: %s" % strands[count]

        print "tsses: %s" % tsses
        print "tsses[count]: %s" % tsses[count]
        print "tss_offsets: %s" % tss_offsets
        print "tss_offsets[count]: %s" % tss_offsets[count]
        lines.append("\t".join((
            str(count),
            avg_pct,
            str(gene_symbols[count]),
            chromosome,
            interval_start,
            interval_end,
            str(tsses[count]),
            str(tss_offsets[count]),
            left_neighbor,
            motif,
            right_neighbor,
            combined_motif,
            cacctg,
            caggtg,
            str(strands[count])
        )))
    return lines




#    db_info = None
#    if options.db_file:
#        db_info = get_db(options)
#    sites = get_sites(options, db_info)
#    (unduplicated_sites, chromosomes) = unduplicate(sites)
#    chromosome_files = read_chrom(chromosomes, options)
#    write_out(unduplicated_sites, chromosome_files, options)
#
#def get_db(options):
#    """create a kvp for symbol : strand"""
#    db_info = {}
#    with open(options.db_file) as db_file:
#        for line in db_file:
#            if line[0] == '#':
#                continue
#
#            line = line.strip()
#            (id, symbol_tss, chrom, start, end, strand) = line.split("\t")
#            (symbol, tss) = symbol_tss.split(":")
#            db_info[symbol] = strand
#    return db_info
#
#def get_sites(options, db_info):
#    """ create sites"""
#    sites = []
#    #open a cpp_result file
#    with open(options.sites_file) as sites_file:
#        #we may have arbitrarily limited results to the first X genes
#        print "limit: %s" % options.limit
#        next(sites_file)
#        for line in sites_file:
#            print "limit: %s, sites: %s\n" % (options.limit, len(sites))
#            if options.limit is not None and len(sites) > options.limit:
#                break
#            #skip the header of a cpp_result
#            if line[0] == "db.desc1":
#                continue
#            #skip the header of a db_file
#            if line[0] == '#':
#                continue
#
#            sites.append(GeneSite(line, db_info))
#
#    #expression rating was used to limit results if limits were desired
#    #now we want to sort for our own purposes as we decide how to deal with overlapping genes
#    sites = sorted(sites, key=lambda site: site.start)
#    sites = sorted(sites, key=lambda site: site.chromosome)
#
#    for site in sites:
#        print site
#
#    return sites
#
#def read_chrom(chromosomes, options):
#    """Read chromosomes in"""
#    print "reading chromosome files."
#
#    chromosome_files = defaultdict(str)
#    #fasta data is separated in chromosome files
#    #if size complexity were a problem, this would be changed to only have one chromosome in memory at a time
#    for chromosome in chromosomes:
#
#        with open(options.chrom_dir + chromosome + ".fa", 'r') as chromosome_file:
#
#            print "read chromosome: %s" % chromosome
#            contents = chromosome_file.read()
#            pos = contents.index('\n')
#            contents = contents[pos+1:]
#            contents = re.sub('[\n\r \t]', '', contents)
#            chromosome_files[chromosome] = contents
#
#
#    print "chromosome files read"
#    return chromosome_files
#
#def unduplicate(sites):
#    """remove duplicate sites and determine chromosomes needed"""
#    unduplicated_sites = []
#    chromosomes = []
#    overlaps = 0
#    pos = 0
#    while pos < len(sites):
#    #    print sites[pos]
#        if sites[pos].chromosome not in chromosomes:
#            chromosomes.append(sites[pos].chromosome)
#
#        pos_start = pos
#        pos_end = pos+1
#        #combine the overlapped sites
#        while pos_end < len(sites) and sites[pos_start].chromosome == sites[pos_end].chromosome and sites[pos_start].end > sites[pos_end].start:
#            overlaps += 1
#            #concatenate to the gene field
#            sites[pos_start].gene += "," + sites[pos_end].gene
#            #concatenate to the tss field
#            sites[pos_start].tss += "," + sites[pos_end].tss
#            sites[pos_start].end = sites[pos_end].end
#            sites[pos_start].strand += "," + sites[pos_end].strand
#            pos_end += 1
#
#        unduplicated_sites.append(sites[pos])
#        pos = pos_end
#
#    print len(sites), "total gene sites."
#    print overlaps, "overlaps."
#    print len(unduplicated_sites), "unduplicated gene sites."
#
#    return (unduplicated_sites, chromosomes)
#
#def write_out(unduplicated_sites, chromosome_files, options):
#    """put it all together now"""
#    #recover the fasta data
#    with open(options.out, 'w') as out:
#        for unduplicated_site in unduplicated_sites:
#            #contents_len = len(chromosome_files[unduplicated_site.chromosome])
#            #print "chromosome length: %d, gene start: %d, gene end: %d" % (contents_len, unduplicated_site.start, unduplicated_site.end)
#
#            snip = chromosome_files[unduplicated_site.chromosome][unduplicated_site.start:unduplicated_site.end]
#
#            if snip == '':
#                print "ERROR: snip was '' for: ", unduplicated_site
#            else:
#                out.write("> %s\n" % unduplicated_site)
#                out.write("%s\n" % snip)
#                out.write("\n")
#
#
if __name__ == "__main__":
    main()
