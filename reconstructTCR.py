import sys
import matcher

if __name__ == '__main__':
	args = matcher.args()

	# sort input file into correct format. Then pass this into matcher.

	matcher.main(args)

	sys.exit()