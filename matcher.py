# matcher.py efficiently finds overlaps between sequences to give fully reconstructed TCR sequences
import sys
import argparse
from time import time

def args():
	parser = argparse.ArgumentParser( description='** script to find overlaps between fragments of TCR sequence and rebuild complete sequences. **')
	parser.add_argument('-f', '--filename', type=str, help='File of sequences to be analysed', required=False)
	parser.add_argument('-mo', '--minoverlap', type=int, help='Minimum length of an overlap', required=False, default=5)
	parser.add_argument('-o', '--outputfile', type=str, help='Write summary to output file', required=False)
	return parser

def getFilename():
		return raw_input("Please enter the name of the file you wish to analyse: ")

def getSequences(file):
	return [line.rstrip() for line in open(file, "r")] # will possibly need to modify as a loop, if file format has more than just sequences

def flatten(l):
  return [item for sublist in l for item in sublist]

def partition(seq, length):
	subseqtree = []
	for i in range(length -1, len(seq)):
		subseqtree.append(seq[len(seq) - 1 - i:])		
	return subseqtree

def search(keyseq,targetseq):
	longest_overlap = ""
	
	for i in range(len(keyseq)):
			if keyseq[i] == targetseq[0:len(keyseq[i])]:
				longest_overlap = keyseq[i]
				break

	return longest_overlap

def searchAll(subseqtree,index,sequences):
	overlaps = []
	for j in range(len(sequences)):
		if index == j:
			continue
		else:
			overlap = search(subseqtree, sequences[j])
		if overlap != "":
			overlaps.append([j,overlap])
	return overlaps

def getFullSeq(seqa, seqb, overlap):
	return seqa + seqb[len(overlap):]

def getSummary(overlaps, sequences, index, targetindex):
	match_index = overlaps[targetindex][0]
	summary =  str(index) + ", " + sequences[index] + "\n"				\
		+ str(match_index) + ", " + sequences[match_index] + "\n"		\
		+ overlaps[targetindex][1] + "\n" + getFullSeq(sequences[index],sequences[match_index],overlaps[targetindex][1]) + "\n\n"
	return summary

def printSummary(overlaps, sequence, index):
	if len(overlaps) == 1:
		for k in range(len(overlaps)):
			print getSummary(overlaps,sequence,index, k)
	elif len(overlaps) > 1:
		print "-------------------------------------------------------------------"
		print "Warning! Multiple overlaps found for same sequence!\n"
		for k in range(len(overlaps)):
			print getSummary(overlaps,sequence,index, k)
		print "-------------------------------------------------------------------"	

def storeSummary(overlaps, sequence, index):
	summary = []
	if len(overlaps) == 1:
		for k in range(len(overlaps)):
			summary.append(getSummary(overlaps,sequence,index, k))
			return summary
	elif len(overlaps) > 1:
		summary.append("-------------------------------------------------------------------\n")
		summary.append("Warning! Multiple overlaps found for same sequence!\n\n")
		for k in range(len(overlaps)):
			summary.append(getSummary(overlaps,sequence,index, k))
		summary.append("-------------------------------------------------------------------\n\n")
		return summary
	else:
		return None

def main(args):
	start_time = time()
	
	if not args.filename: args.filename = getFilename()  # possibly build function to validate file

	sequences = getSequences(args.filename)
	stored_summaries = []
	
	for i in range(len(sequences)):
		subseqtree = partition(sequences[i],args.minoverlap)
		overlaps = searchAll(subseqtree,i,sequences)
		if args.outputfile:
			summary = storeSummary(overlaps,sequences,i)
			if summary:
				stored_summaries.append(summary)
		else:
			printSummary(overlaps,sequences,i)

	if args.outputfile:
		file = open(args.outputfile,'w')
		file.write("".join(flatten(stored_summaries)))
		file.close()
		print "Output written to",args.outputfile

	print "matching complete after ", (time() - start_time)


if __name__ == '__main__':
	parser = args()
	args = parser.parse_args()
	main(args)
	sys.exit()