import sys
import re

# line = "##contig=<ID=scaffold941|size354645|quiver,length=355866>"
# line = "scaffold45|size1037999|quiver   167436  .       C       A       5.76485 .       DP=13;VDB=0.78;SGB=-2.48712;RPB=0.85;MQB=0.1;MQSB=0.594978;BQB=0.9;MQ0F=0;ICB=0.5;HOB=0.5;AC=2;AN=4;DP4=4,6,0,2;MQ=55   GT:PL:AD  0/1:31,0,31:1,1  0/1:7,0,219:9,1"

infoCols = 8
startIndCols = infoCols + 1

freq = open("refFreqs.txt", "w")
count = open("totCounts.txt", "w")

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

with open('alfalfaAlignCADL/chmf.ploidy.vcf') as f:
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


with open('alfalfaAlignCADL/chmf.ploidy.vcf') as f:
# with open('test.vcf') as f:
	for line in f:
		line = line.split('#', 1)[0]
		line = line.rstrip()
		# if not li.startswith("#"):
		if line:
			cntaf = map(getAf, line.split()[startIndCols:])
			# print(af)
			cnt = [el[0] for el in cntaf]
			af = [el[1] for el in cntaf]
			freq.write(" ".join(af) + "\n")
			count.write(" ".join(cnt) + "\n")
	freq.write('\n')
	count.write('\n')
		# else:
			# print("this line is useless!", line)

freq.close()
count.close()

