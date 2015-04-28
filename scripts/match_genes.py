import sys
import re
from collections import defaultdict

class GeneSite:
    def __init__(self, line):
#        print line
        words = re.split("[ \t]", line)
#        print words
        (symbol, tss) = words[1].split(":")
        self.gene = symbol
        self.chromosome = words[2]
        self.start = int(words[3])
        self.end = int(words[4])

#        words = re.split("[ \t]", line)
#        self.gene = words[1]
#        self.chromosome = words[2]
#        self.start = int(words[3])
#        self.end = int(words[4])

    def __repr__(self):
        return "<%s %s %s %s>" % (self.gene, self.chromosome, self.start, self.end)

    def __str__(self):
        return "'%s %s %s %s'" % (self.gene, self.chromosome, self.start, self.end)


sites = []
with open(sys.argv[1]) as sites_file:
    print sys.argv
    if len(sys.argv) > 4:
        limit = int(sys.argv[4])
    else:
        limit = None
    print "limit: %s" % limit
    next(sites_file)
    for line in sites_file:
        print "limit: %s, sites: %s\n" % (limit, len(sites))
        if limit is not None and len(sites) > limit:
            break
        if line[0] == "db.desc1":
            continue
        if line[0] == '#':
            continue
        site = GeneSite(line)
        if site.start > site.end:
            print "need to swap start/end for:", site
            tmp = site.start
            site.start = site.end
            site.end = tmp

        sites.append(site)

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
#        print "overlap detected between:";
#        print "\t", sites[x_start];
#        print "\t", sites[x_end];
        overlaps += 1
        sites[x_start].gene += ", " + sites[x_end].gene
        sites[x_start].end = sites[x_end].end
        x_end += 1

#    print "appending:", sites[x]
    unduplicated_sites.append(sites[x])
    x = x_end

print len(sites), "total gene sites."
print overlaps, "overlaps."
print len(unduplicated_sites), "unduplicated gene sites."

print "reading chromosome files."

chromosome_files = defaultdict(str)
for chromosome in chromosomes:
    with open(sys.argv[2] + chromosome + ".fa", 'r') as chromosome_file:

        contents = chromosome_file.read()
        pos = contents.index('\n')
        contents = contents[pos+1:]
        contents = re.sub('[\n\r \t]', '', contents)
#        print "contents[0:1000]: '" + contents[0:1000] + "'"

        chromosome_files[chromosome] = contents

        print "read chromosome", chromosome

print "read chromosome files.";

with open(sys.argv[3], 'w') as out:
    for unduplicated_site in unduplicated_sites:
        contents_len = len(chromosome_files[unduplicated_site.chromosome])
        #print "chromosome length: %d, gene start: %d, gene end: %d" % (contents_len, unduplicated_site.start, unduplicated_site.end)

        snip = chromosome_files[unduplicated_site.chromosome][unduplicated_site.start:unduplicated_site.end]

        if snip == '':
             print "ERROR: snip was '' for: ", unduplicated_site
        else:
             out.write("> %s\n" % unduplicated_site)
             out.write("%s\n" % snip)
             out.write("\n")
#       print ">", unduplicated_site
#       print snip
#       print ""
