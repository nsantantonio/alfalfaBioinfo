import sys
import re
infile = sys.argv[1]

# lookup = 'CHROM'
with open(infile) as f:
    for i, line in enumerate(f, 1):
        if re.match("\\#CHROM.*", line):
            print 'First line of scores in vcf:', i + 1
    		break
