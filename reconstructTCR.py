import os
import sys
import matcher
import argparse
from IPython import embed

def addargs(parser):
	parser.add_argument('-nc', '--nocollapse', action='store_true', help='Prevents sequences being with same tag being collapsed into longest possible sequence', required=False, default=False)
	parser.add_argument('-nr', '--noreconstruct', action='store_true', help='Prevents reconstruct of collapsed sequences into full TCRs', required=False, default=False)
	parser.add_argument('-s', '--separatedir', type=str, help='Name of directory in which to store tags separated into unique files', required=False)	
	parser.add_argument('-c', '--chains', type=str, help='Specify chains if known', nargs='*', required=False, default=['a','b','c','d'])
	return parser

def argchecker():

	if args.nocollapse and not args.noreconstruct:
		print "Warning: It is inadivsable to reconstruct TCRs if input sequences have not been collapsed.\n" \
			  "Only use the '-nc' argument if input sequences have already been collapsed.\n"

def separateIntoFiles():
	chain_names = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}
	parent = args.separatedir

	if not os.path.exists(parent):
		os.makedirs(parent)
	else:
		overwrite = raw_input("Warning! Directory '"+parent+"' already exists. Do you wish to overwrite? (Y/N): ")
		if overwrite not in ["Y","y","yes"]:
			print "Directory '"+parent+"' was not overwritten. Program will now exit."
			sys.exit()

	for tag in separated_seqs.keys():
		
		for chain in separated_seqs[tag].keys():
			if separated_seqs[tag][chain]:
				if chain in chain_names:
					dir_name = parent+"/"+chain_names[chain]
				else:
					dir_name = chain
				if not os.path.exists(dir_name):
					os.makedirs(dir_name)

				for gene in separated_seqs[tag][chain]:
					tag_file = dir_name+"/"+tag+gene
					file = open(tag_file,'w')
					for dcr in separated_seqs[tag][chain][gene]:
						file.write(dcr)
					file.close()

	print "Separated files can be found in", parent


def specifiedChain(chain):
	if chain in args.chains:
		return True

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
		if not specifiedChain(chain):
			continue

		elif dcr[1] != "n/a" and dcr[2] == "n/a":
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

def getIndex(dictionary, search_value):
	return [key for key, val in dictionary.iteritems() if val == search_value][0]

def getTargetIndicies(overlaps, chain, tag, dictionary, searchseqs):
	target_indicies = []
	for i in range(len(overlaps)):
		target_indicies.append(chain+tag+getIndex(dictionary,searchseqs[overlaps[i][0]]))
	return target_indicies

def annotateSummary(match, search_index, target_indicies):

	if len(match) == 1:
		match_range = range(1)
		ishift = 0
	else:
		match_range = range(2,len(match)-1)
		ishift = 2

	for i in match_range:
		splitmatch = match[i].split("\n")
		old_search_index = splitmatch[0][0]
		old_target_index = splitmatch[1][0]
		splitmatch[0] = splitmatch[0].replace(old_search_index,search_index)
		splitmatch[1] = splitmatch[1].replace(old_target_index,target_indicies[i-ishift])
		newrecord = "\n".join(splitmatch)
		match[i] = newrecord
	
	return match

def findMatches():
	summary = []

	for tag_type in [['v','j'],['j','v']]:
		
		for chain in longest_seqs[tag_type[0]].keys():
			
			for key in longest_seqs[tag_type[0]][chain].keys():

				seq = longest_seqs[tag_type[0]][chain][key]
				if chain in longest_seqs[tag_type[1]].keys():
					searchseqs = [seq] + longest_seqs[tag_type[1]][chain].values()
					subseqtree = matcher.partition(seq,args.minoverlap)
					overlaps = matcher.searchAll(subseqtree,0,searchseqs)
					search_index = chain+tag_type[0]+key
					target_indicies = getTargetIndicies(overlaps, chain, tag_type[1], longest_seqs[tag_type[1]][chain], searchseqs)
					match = matcher.storeSummary(overlaps,searchseqs,0)
				
				if match:
					summary.append(annotateSummary(match,search_index,target_indicies))
	return summary

def printDuplicates(summary):
	for i in range(len(summary)):
		if len(summary[i]) > 1:
			for j in range(len(summary[i])):
				print summary[i][j]

def writeSummary(summary):
	summary_for_file = []
	
	for i in range(len(summary)):
		if len(summary[i]) == 1:
			summary_for_file.append(summary[i])
		else:
			summary_for_file.append(summary[i][2:len(summary[i])-1])
	
		file = open(args.outputfile,'w')
		file.write("".join(matcher.flatten(summary_for_file)))
		file.close()

	printDuplicates(summary)
	if summary:
		print "Output written to",args.outputfile
	else:
		print "No matches found!"
	
def printSummary(summary):
	print "---------------------------------------------------------------------------"
	print "                                 Matches                                   "
	print "---------------------------------------------------------------------------"
	for i in range(len(summary)):
		for j in range(len(summary[i])):
			print summary[i][j]
	if not summary:
		print "No matches found!"
	
if __name__ == '__main__':
	 
	parser = matcher.args()
	args = addargs(parser).parse_args()
	argchecker()

	chains = args.chains

	v_tags, j_tags, double_tags = separate(chains)

	separated_seqs = {'v': v_tags, 'j': j_tags, 'vj': double_tags}

	if args.separatedir:
		separateIntoFiles()

	if not args.nocollapse:
		longest_seqs = {}
		for tag_type in separated_seqs.keys():
			longest_seqs[tag_type] = collapseTag(separated_seqs[tag_type], tag_type)
	
	if not args.noreconstruct:
		summary = findMatches()

		if args.outputfile:
			writeSummary(summary)
		else:
			printSummary(summary)

	sys.exit()