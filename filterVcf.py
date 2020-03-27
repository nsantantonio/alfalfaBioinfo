import time
startTime = time.time()

import sys
import re
import argparse



parser = argparse.ArgumentParser(description='Calculate allele frequencies from vcf and filter based on min, mean and total read counts per site')
parser.add_argument("-f", '--infile', type=str, help='vcf file name')
parser.add_argument("-o", '--outdir', type=str, default = "", help='output directory and file prefix')
parser.add_argument('--minCnt', type=int, default=0,
                    help='Minimum read count for all genotyped lines. Any sites with at least one genotype with a read count below this integer will be excluded')
parser.add_argument('--maxCnt', type=int, default=200,
                    help='Maximum read count for all genotyped lines. Any sites with at least one genotype with a read count above this integer will be excluded')
parser.add_argument('--minMean', type=int, default=20,
                    help='Minimum mean read count. Any sites with a mean read count below this integer will be excluded')
parser.add_argument('--maxMean', type=int, default=200,
                    help='Maximum mean read count. Any sites with a mean read count above this integer will be excluded')
parser.add_argument('--minTot', type=int, default=0,
                    help='Minimum total read count. Any sites with a total read count below this integer will be excluded')
parser.add_argument('--maxTot', type=int, default=10000,
                    help='Maximum total read count. Any sites with a total read count above this integer will be excluded')

args = parser.parse_args()

infile = args.infile
outdir = args.outdir
minCnt = args.minCnt
maxCnt = args.maxCnt
minMean = args.minMean
maxMean = args.maxMean
minTot = args.minTot
maxTot = args.maxTot

# line = "Super-Scaffold_1        357018  .       C       T       180     .       DP=38;VDB=0.0278409;SGB=2.66165;RPB=0.67934;MQB=0.990262;MQSB=0.93828;BQB=0.714832;MQ0F=0;ICB=1;HOB=0.0798611;AC=5;AN=24;DP4=11,18,3,4;MQ=59       GT:PL:AD        ./.:0,0,0:0,0   0/0:0,3,37:1,0  0/0:0,3,37:1,0  ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   0/1:51,0,65:2,1 ./.:0,0,0:0,0   ./.:0,0,0:0,0   0/0:0,3,25:1,0  0/0:0,3,37:1,0  1/1:105,9,0:0,3    ./.:0,0,0:0,0   0/1:67,6,0:0,2  ./.:0,0,0:0,0   0/0:0,3,37:1,0  ./.:0,0,0:0,0   ./.:0,0,0:0,0   0/0:0,9,92:3,0  ./.:0,0,0:0,0   0/0:0,39,255:13,0       0/1:19,0,134:5,1        ./.:0,0,0:0,0   0/0:0,3,37:1,0"
# infile = "test/test.vcf"
# outdir = ""
# minCnt = 0
# maxCnt = 200
# minMean = 1
# maxMean = 200
# minTot = 10
# maxTot = 10000

infoCols = 8
startIndCols = infoCols + 1
infokeep = 6

suffix = "_cnt" + str(minCnt) + "-" + str(maxCnt) + "_mean" + str(minMean) + "-" + str(maxMean) + "_tot" + str(minTot) + "-" + str(maxTot)

info = open(outdir + "namePos" + suffix + ".txt", "w")
freq = open(outdir + "refFreqs" + suffix + ".txt", "w")
count = open(outdir + "totCounts" + suffix + ".txt", "w")

def mean(x):
	mu = sum(x) / float(len(x))
	return mu

def getAf(x, filter = True):
	y = x.rsplit(":",1)[1]
	a = [int(z) for z in y.split(',')]
	suma = sum(a)
	if suma == 0:
		r = [0, 0]
	else:
		r = [suma, a[0] / float(suma)]
	return r


with open(infile) as f:
# with open('test.vcf') as f:
	for line in f:
		if re.match("\\#CHROM.*", line):
			infonames = line.split()[0:infokeep]
			names = [re.sub('.*\\/|\\.bam', '', i) for i in line.split()[startIndCols:]]
		line = line.split('#', 1)[0]
		line = line.rstrip()
		if line:
			n = len(line.split()[startIndCols:])
			break
# print(infonames)
infonames = 
info.write(" ".join(infonames) + "\n")
freq.write(" ".join(names) + "\n")
count.write(" ".join(names) + "\n")

with open(infile) as f:
	for line in f:
		line = line.split('#', 1)[0]
		line = line.rstrip()
		if line:
			scafPos = line.split()[0:infokeep]
			# x = line.split()[startIndCols+18]
			cntaf = map(getAf, line.split()[startIndCols:])
			# check if python 2 or 3
			if (sys.version_info > (3, 0)):
				cntaf = list(cntaf)
			cnt = [i[0] for i in cntaf]
			af = [i[1] for i in cntaf]
			minl = min(cnt)
			maxl = max(cnt)
			totCnt = sum(cnt)
			meanCnt = mean(cnt)
			if minl >= minCnt and maxl <= maxCnt and meanCnt >= minMean and meanCnt <= maxMean and totCnt >= minTot and totCnt <= maxTot :
				info.write(" ".join(scafPos) + "\n")
				freq.write(" ".join([str(i) for i in af]) + "\n")
				count.write(" ".join([str(i) for i in cnt]) + "\n")
	info.write('\n')
	freq.write('\n')
	count.write('\n')
info.close()
freq.close()
count.close()

secs = time.time() - startTime
if secs > 3600:
	t = round(secs / 3600, 1)
	print("Run time: %s hours" % (t))
elif secs > 60:
	t = round(secs / 60, 1)
	print("Run time: %s minutes" % (t))
else:
	t = round(secs, 3)
	print("Run time: %s seconds" % (t))
