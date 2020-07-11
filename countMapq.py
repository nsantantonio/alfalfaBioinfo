import sys
import collections

with open(sys.argv[0]) as infile:
    counts = collections.Counter(l.strip() for l in infile)
for line, count in counts.most_common():
    print line, count