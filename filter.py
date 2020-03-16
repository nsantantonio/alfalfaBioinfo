import sys
import re

# line = "##contig=<ID=scaffold941|size354645|quiver,length=355866>"
# line = "scaffold45|size1037999|quiver   167436  .       C       A       5.76485 .       DP=13;VDB=0.78;SGB=-2.48712;RPB=0.85;MQB=0.1;MQSB=0.594978;BQB=0.9;MQ0F=0;ICB=0.5;HOB=0.5;AC=2;AN=4;DP4=4,6,0,2;MQ=55   GT:PL:AD  0/1:31,0,31:1,1  0/1:7,0,219:9,1"

# line = "Super-Scaffold_1        357018  .       C       T       180     .       DP=38;VDB=0.0278409;SGB=2.66165;RPB=0.67934;MQB=0.990262;MQSB=0.93828;BQB=0.714832;MQ0F=0;ICB=1;HOB=0.0798611;AC=5;AN=24;DP4=11,18,3,4;MQ=59       GT:PL:AD        ./.:0,0,0:0,0   0/0:0,3,37:1,0  0/0:0,3,37:1,0  ./.:0,0,0:0,0   ./.:0,0,0:0,0   ./.:0,0,0:0,0   0/1:51,0,65:2,1 ./.:0,0,0:0,0   ./.:0,0,0:0,0   0/0:0,3,25:1,0  0/0:0,3,37:1,0  1/1:105,9,0:0,3    ./.:0,0,0:0,0   0/1:67,6,0:0,2  ./.:0,0,0:0,0   0/0:0,3,37:1,0  ./.:0,0,0:0,0   ./.:0,0,0:0,0   0/0:0,9,92:3,0  ./.:0,0,0:0,0   0/0:0,39,255:13,0       0/1:19,0,134:5,1        ./.:0,0,0:0,0   0/0:0,3,37:1,0"
# infile = "alfalfaVariantCalls/test.vcf"
# outdir = ""

infile = sys.argv[1]
outdir = sys.argv[2]


fprefix = "test/test_"
# outdir = sys.argv[2]

namePos
refFreqs
totCounts


def mean(x):
	mu = sum(x) / len(x)
	return mu

# from itertools import izip

with open(fprefix + "namePos") as counts, open("refFreqs") as freqs: 
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

