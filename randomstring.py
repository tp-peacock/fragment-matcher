# generate random matching and non matching strings
import sys
import random, string
import argparse
from time import time
from IPython import embed

def args():
	parser = argparse.ArgumentParser( description='** generates file of random strings, pairs of matching strings, matching strings with deleted characters **')
	parser.add_argument('-t', '--total', type=int, help='Total number of strings', required=False, default=100)
	parser.add_argument('-l', '--maxlength', type=int, help='Max length of string', required=False, default=50)
	parser.add_argument('-o', '--maxoverlap', type=int, help='Max length of string overlap', required=False, default=12)
	parser.add_argument('-p', '--matchchance', type=int, help='Percentage chance of two matching strings', required=False, default=10)
	parser.add_argument('-f', '--outputfile', type=str, help='Output filename (required)', required=True)
	return parser.parse_args()

def randomString(length):
   return ''.join(random.choice(string.lowercase) for i in range(length))

def reverse(string):
	return string[::-1]

if __name__ == '__main__':
	start_time = time()
	args = args()

	strings = []

	string_count = 0
	paired_strings_count = 0

	while string_count < args.total:
		randstring = randomString(args.maxlength - args.maxoverlap)
		chance = random.randint(1,100)
		if chance <= args.matchchance and string_count <= args.total - 2:
			overlap = randomString(args.maxoverlap)
			pair_a = randstring+overlap
			pair_b = overlap+(randomString(args.maxlength - args.maxoverlap))
			strings.append(pair_a)
			strings.append(pair_b)
			string_count += 2
			paired_strings_count +=1
		else:
			strings.append(randstring+randomString(args.maxoverlap))
			string_count += 1
	random.shuffle(strings)

	file = open(args.outputfile,"w")
	for i in range(len(strings)):
		file.write(strings[i]+"\n")
	file.close()

	# SUMMARY
	print "Strings written to", args.outputfile
	print "Total number of strings is ", string_count
	print "Total number of paired strings is ", paired_strings_count
	print "Generation completed after ", (time() - start_time)

	sys.exit() 	 	