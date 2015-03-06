#!/usr/bin/python
import MySQLdb
import re
from optparse import OptionParser

#usage = "usage: %prog [options] arg"
usage = ""
parser = OptionParser(usage)
parser.add_option("--step", dest="step", help="Which step should be compared?", metavar="INT")
parser.add_option("--sample_name", dest="sample_name", help="What sample should be used?", metavar="STRING")
parser.add_option("--passwd", dest="passwd", help="What passwd to use", metavar="STRING")

(options, args) = parser.parse_args()
if not options.step or not options.sample_name:
	parser.print_help()
	parser.error("step or sample_name not given")


db = MySQLdb.connect(
		host="localhost", 
		user="tdesell",
		passwd=options.passwd,
		db="wildlife")

cur = db.cursor()
cur.execute(
	"\n".join((
		"SELECT name, sequences_filename, command_line_options, gibbs_walk.id, gibbs_walk.current_steps",
		"FROM gibbs_sampler JOIN gibbs_walk ON sampler_id = gibbs_sampler.id",
		"WHERE name like %s"
		)),
	options.sample_name
)

fields = map(lambda x:x[0], cur.description)
result = [dict(zip(fields, row)) for row in cur.fetchall()]
if len(result) == 0:
	parser.print_help()
	parser.error("sample_name not found")

header_r = re.compile("^[FR]|$")
body_r = re.compile("^\s+(\d+),\s*(\d+)\s+(\d+)\s+([a-z]*)([A-Z]+)([a-z]*)\s+(\d+)\s+([\d\.]+)\s+>\s+'([^']+)'$")
#23, 1      31 tgtgcTCTCCTaaacc      36 0.1702 > 'PNP chr14 20937037 20938037'



motifs = {}

def process_line(line):
	
	if header_r.match(line):
		return
	else:
		body_m = body_r.match(line)
		if body_m is not None:
			(seq, motif_num, start, lneighbor, motif, rneighbor, end, pct, id) = body_m.groups()
			if motif not in motifs:
				motifs[motif] = {}
			if id not in motifs[motif]:
				motifs[motif][id] = {"lneighbor":lneighbor, "rneighbor":rneighbor, "pct":0.0}
			motifs[motif][id]["pct"] += float(pct)
		else:
			print "no match? %s" % line

#motif_files = ["./walk_200010_steps_30000.motifs"]
#for motif_file in motif_files:
for walk in result:

	motif_file = "/data/dna_at_home/%s/walk_%s_steps_%s.motifs" % (walk["name"], walk["id"], options.step)
#	print "new : %s\n" % new_args
	with open(motif_file, "r") as in_file: 
		try:
			gen = (process_line(line) for line in in_file)

			while 1:
				next(gen)
			
		except StopIteration:
			print "file %s done" % motif_file
			in_file.close()
		
print "motifs: %s" % motifs
