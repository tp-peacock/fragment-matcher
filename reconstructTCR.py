import sys
import matcher
import argparse
from IPython import embed

def addargs(parser):
	parser.add_argument('-nc', '--nocollapse', action='store_true', help='Prevents sequences being with same tag being collapsed into longest possible sequence', required=False, default=False)
	parser.add_argument('-nr', '--noreconstruct', action='store_true', help='Prevents reconstruct of collapsed sequences into full TCRs', required=False, default=False)
	parser.add_argument('-s', '--separate', action='store_true', help='Separate tags into separate files', required=False, default=False)	
	return parser

def argchecker():

	if args.nocollapse and not args.noreconstruct:
		print "Warning: It is inadivsable to reconstruct TCRs if input sequences have not been collapsed.\n" \
			  "Only use the '-nc' argument if input sequences have already been collapsed.\n"

def separate(chains):
	f = open(args.filename, 'r')
	v_tags = {}
	j_tags = {}
	double_tags = {}

	for i in range(len(chains)):
		v_tags[chains[i]] = {}
		j_tags[chains[i]] = {}
		double_tags[chains[i]] = {}

	for line in f:
		dcr = line.split(", ")

		chain = dcr[0]
	
		if dcr[1] != "n/a" and dcr[2] == "n/a":
			if dcr[1] in v_tags[chain]:
				v_tags[chain][dcr[1]].append(line)
			else:
				v_tags[chain][dcr[1]] = []
				v_tags[chain][dcr[1]].append(line)

		elif dcr[1] == "n/a" and dcr[2] != "n/a":
			if dcr[2] in j_tags[chain]:
				j_tags[chain][dcr[2]].append(line)
			else:
				j_tags[chain][dcr[2]] = []
				j_tags[chain][dcr[2]].append(line)

		elif dcr[1] != "n/a" and dcr[2] != "n/a":
			print "--------------------------------------------------------------------------------------------"
			print "double tag found in: "+args.filename
			print dcr
			print "--------------------------------------------------------------------------------------------"
			if ("v"+dcr[1]+"_j"+dcr[2]) in double_tags[chain]:
				double_tags[chain]["v"+dcr[1]+"_j"+dcr[2]].append(line)
			else:
				double_tags[chain]["v"+dcr[1]+"_j"+dcr[2]] = []
				double_tags[chain]["v"+dcr[1]+"_j"+dcr[2]].append(line)
	
	f.close()
	return v_tags, j_tags, double_tags

def mostCommon(lst):
    return max(set(lst), key=lst.count)

def collapseTag(dcrs, tag_type):
	longest_seqs = {}
	for chain in dcrs.keys():	
		if dcrs[chain] != {}:
			longest_seqs[chain] = collapseChain(dcrs[chain], tag_type)
	return longest_seqs

def collapseChain(tags, tag_type):
	longest_seqs = {}
	for tag in tags.keys():
		longest_seqs[tag] = collapse(tags[tag],tag_type)
	return longest_seqs

def collapse(dcrs, tag_type):
	positions = {}
	for i in range(len(dcrs)):
		seq = dcrs[i].split(", ")[4]
		if tag_type == "j":
			seq = seq[::-1]
		for base in range(len(seq)):
			if base in positions:
				positions[base].append(seq[base])
			else:
				positions[base]=[]
				positions[base].append(seq[base])
	longest_seq = buildLongestSeq(positions)
	if tag_type == "j":
		longest_seq = longest_seq[::-1]
	return longest_seq

def buildLongestSeq(positions):
	seq_bases = []
	for base in positions:
		seq_bases.append(mostCommon(positions[base]))
	longest_seq = ("").join(seq_bases)
	return longest_seq


if __name__ == '__main__':
	 
	parser = matcher.args()
	args = addargs(parser).parse_args()
	argchecker()

	chains= ['a','b','c','d']
	v_tags, j_tags, double_tags = separate(chains)
	separated_seqs = {'v': v_tags, 'j': j_tags, 'vj': double_tags}

	if args.separate:
		print "Unfinished: program does not yet fully support separation\n" #separate tags here

	if not args.nocollapse:
		longest_seqs = {}
		for tag_type in separated_seqs.keys():
			longest_seqs[tag_type] = collapseTag(separated_seqs[tag_type], tag_type)

	print longest_seqs
	embed()

	# sort input file into correct format. Then pass this into matcher.

	# matcher.main(args)

	sys.exit()