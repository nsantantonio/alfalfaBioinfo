import sys
import re

import argparse

parser = argparse.ArgumentParser(description='filter vcf based on min, mean and total read counts per site')

parser.add_argument("-f", '--prefix', type=str, default = "", help='vcf file name prefix')

parser.add_argument('--minCnt', type=int, default=0,
                    help='Minimum read count for all genotyped lines. Any sites with at least one genotype with a read count below this integer will be excluded')
parser.add_argument('--maxCnt', type=int, default=0,
                    help='Maximum read count for all genotyped lines. Any sites with at least one genotype with a read count above this integer will be excluded')


parser.add_argument('--minMn', type=int, default=20,
                    help='Minimum mean read count. Any sites with a mean read count below this integer will be excluded')
parser.add_argument('--maxMn', type=int, default=20,
                    help='Maximum mean read count. Any sites with a mean read count above this integer will be excluded')


parser.add_argument('--minTot', type=int, default=100,
                    help='Minimum total read count. Any sites with a total read count below this integer will be excluded')
parser.add_argument('--maxTot', type=int, default=100,
                    help='Maximum total read count. Any sites with a total read count above this integer will be excluded')

args = parser.parse_args()

print(args.fprefix + "_cnt" + str(args.minCnt) + "-" + str(args.maxCnt) + "_" + "mean" + str(args.minMn) + "-" + str(args.maxMn) + "_" + "tot" + str(args.minTot) + "-" + str(args.maxTot))


sys.exit("Gave up trying to read two files simultaneously. use filterVcf.py instead to filter the vcf file directly. Not ideal-")
# line = "Super-Scaffold_1        1802    .       G       C       4.06312 .       DP=2;SGB=-0.0106344;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=31  GT:PL:AD        ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0      ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0      1/1:25,3,0:0,1  ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0"
# infile = "test/test.vcf"
# outdir = ""

# infile = 
# outdir = sys.argv[2]
# minCnt = 0
# maxCnt = 0
# minMn = 20
# maxMn = 20
# minTot = 100
# maxTot = 100
# fprefix = "test/test_"

def mean(x):
	mu = sum(x) / len(x)
	return mu

# from itertools import izip

with open(fprefix + "totCounts.txt") as counts, open("refFreqs.txt") as freqs: 
    for x, y in zip(counts, freqs):
        x = x.strip()
        y = y.strip()
        print("{0}\t{1}".format(x, y))


from itertools import izip

with open("textfile1") as textfile1, open("textfile2") as textfile2: 
    for x, y in izip(textfile1, textfile2):
        x = x.strip()
        y = y.strip()
        print("{0}\t{1}".format(x, y))

info = open(outdir + "namePos.txt", "w")
freq = open(outdir + "refFreqs.txt", "w")
count = open(outdir + "totCounts.txt", "w")

def getAf(x):
	y = x.rsplit(":",1)[1]
	a = [int(z) for z in y.split(',')]
	suma = sum(a)
	if suma == 0:
		r = ["0", "0"]
	else:
		r = [str(suma), str(a[0] / float(suma))]
		# r = str(a[0] / float(suma))
	return r

with open(infile) as f:
# with open('test.vcf') as f:
	for line in f:
		if re.match("\\#CHROM.*", line):
			names = [re.sub('.*\\/|\\.bam', '', i) for i in line.split()[startIndCols:]]
		line = line.split('#', 1)[0]
		line = line.rstrip()
		if line:
			n = len(line.split()[startIndCols:])
			break

freq.write(" ".join(names) + "\n")
count.write(" ".join(names) + "\n")






# result.write("Count refFreq\n")

# freqs = [name + "_freq" for name in names]
# counts = [name + "_cnt" for name in names]


with open(infile) as f:
# with open('test.vcf') as f:
	for line in f:
		line = line.split('#', 1)[0]
		line = line.rstrip()
		# if not li.startswith("#"):
		if line:
			scafPos = line.split()[0:5]
			cntaf = map(getAf, line.split()[startIndCols:])
			# print(af)
			cnt = [el[0] for el in cntaf]
			af = [el[1] for el in cntaf]
			info.write(" ".join(scafPos) + "\n")
			freq.write(" ".join(af) + "\n")
			count.write(" ".join(cnt) + "\n")
	freq.write('\n')
	count.write('\n')
		# else:
			# print("this line is useless!", line)

info.close()
freq.close()
count.close()

