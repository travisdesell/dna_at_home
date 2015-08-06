import sys
import re
from collections import defaultdict

#every line from the cpp_result input file is a new gene site
class GeneSite:
    def __init__(self, line):
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
        #determine a strand based on the direction of the region
        self.strand = "+"
        if self.start > self.end:
            #adjust for the direction of the region, gibbs sampler wants everything in the same direction
            #it will guess at directionality
            print "need to swap start/end for:", site
            tmp = self.start
            self.start = self.end
            self.end = tmp
            self.strand = "-"

    def __repr__(self):
        return "<%s %s %s %s %s %s>" % (self.gene, self.chromosome, self.start, self.end, self.tss, self.strand)

    def __str__(self):
        return "'%s %s %s %s %s %s'" % (self.gene, self.chromosome, self.start, self.end, self.tss, self.strand)


sites = []
#open a cpp_result file
with open(sys.argv[1]) as sites_file:
    print sys.argv
    if len(sys.argv) > 4:
        limit = int(sys.argv[4])
    else:
        limit = None
    #we may have arbitrarily limited results to the first X genes
    print "limit: %s" % limit
    next(sites_file)
    for line in sites_file:
        print "limit: %s, sites: %s\n" % (limit, len(sites))
        if limit is not None and len(sites) > limit:
            break
        #skip the header of a cpp_result
        if line[0] == "db.desc1":
            continue
        #skip the header of a db_file
        if line[0] == '#':
            continue
        
        sites.append(GeneSite(line))

#expression rating was used to limit results if limits were desired
#now we want to sort for our own purposes as we decide how to deal with overlapping genes
sites = sorted(sites, key=lambda site: site.start)
sites = sorted(sites, key=lambda site: site.chromosome)

for site in sites:
    print site


unduplicated_sites = []
chromosomes = []
overlaps = 0
x = 0
while x < len(sites):
#    print sites[x]
    if sites[x].chromosome not in chromosomes:
        chromosomes.append(sites[x].chromosome)

    x_start = x
    x_end = x+1
    #combine the overlapped sites
    while x_end < len(sites) and sites[x_start].chromosome == sites[x_end].chromosome and sites[x_start].end > sites[x_end].start:
        overlaps += 1
        #concatenate to the gene field
        sites[x_start].gene += "," + sites[x_end].gene
        #concatenate to the tss field
        sites[x_start].tss += "," + sites[x_end].tss
        sites[x_start].end = sites[x_end].end
        x_end += 1

    unduplicated_sites.append(sites[x])
    x = x_end

print len(sites), "total gene sites."
print overlaps, "overlaps."
print len(unduplicated_sites), "unduplicated gene sites."

print "reading chromosome files."

chromosome_files = defaultdict(str)
#fasta data is separated in chromosome files
#if size complexity were a problem, this would be changed to only have one chromosome in memory at a time
for chromosome in chromosomes:

    with open(sys.argv[2] + chromosome + ".fa", 'r') as chromosome_file:

        contents = chromosome_file.read()
        pos = contents.index('\n')
        contents = contents[pos+1:]
        contents = re.sub('[\n\r \t]', '', contents)
        chromosome_files[chromosome] = contents

        print "read chromosome", chromosome

print "read chromosome files.";

#recover the fasta data
with open(sys.argv[3], 'w') as out:
    for unduplicated_site in unduplicated_sites:
        contents_len = len(chromosome_files[unduplicated_site.chromosome])
        #print "chromosome length: %d, gene start: %d, gene end: %d" % (contents_len, unduplicated_site.start, unduplicated_site.end)

        snip = chromosome_files[unduplicated_site.chromosome][unduplicated_site.start:unduplicated_site.end]

        if snip == '':
             print "ERROR: snip was '' for: ", unduplicated_site
        else:
            #fix indent
             out.write("> %s\n" % unduplicated_site)
             out.write("%s\n" % snip)
             out.write("\n")
#       print ">", unduplicated_site
#       print snip
#       print ""
