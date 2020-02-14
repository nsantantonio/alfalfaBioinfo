nC=40
mapq=20
genome="/local/rdata/dbgap_robbins/alfalfa3a.2/3a.2-bionano+unscaffolded-assembly.fa"
# using just the first 978 scaffolds, as they appear to be the main genome sequence. Wrote Andrew Farmer email to confirm
# genome="/local/rdata//local/rdata/dbgap_robbins/alfalfaMainGenome/3a.2-bionano+unscaffolded-assembly_1to978.fa"
workDir=alfalfaWorkdir
testDir=alfalfaShortFastqTest
seqDir=alfalfaNovaSeq
# bamDir=alfalfaBams
bamDir=testBams
ls ${seqDir}/*.fastq.gz > ${workDir}/seqFiles.txt
sed 's/_L001_.*\|_L002_.*\|.*\///g' ${workDir}/seqFiles.txt | uniq > ${workDir}/gidWells.txt 
sed 's/-Rep._.*//g' ${workDir}/gidWells.txt | sed 's/.*-//g' | uniq > ${workDir}/gid.txt 
sed 's/.*alfalfaDNAext-//g' ${workDir}/gidWells.txt | sed 's/-.*//g' | uniq > ${workDir}/wells.txt 
cat ${workDir}/seqFiles.txt
cat ${workDir}/gidWells.txt
cat ${workDir}/gid.txt
cat ${workDir}/wells.txt

# head -2 ${workDir}/gidWells.txt > ${workDir}/test.txt

while read f; do
	newf="${f##*/}"
	zcat $f | head -1000 | gzip > ${testDir}/${newf}
done < ${workDir}/seqFiles.txt

while read g; do
	bwa mem -t ${nC} ${genome} ${testDir}/${g}_L001_R1_001.fastq.gz ${testDir}/${g}_L001_R2_001.fastq.gz > ${bamDir}/${g}_L001_aln.sam
	bwa mem -t ${nC} ${genome} ${testDir}/${g}_L002_R1_001.fastq.gz ${testDir}/${g}_L002_R2_001.fastq.gz > ${bamDir}/${g}_L002_aln.sam
	samtools view -b -q ${mapq} -@${nC} ${bamDir}/${g}_L001_aln.sam > ${bamDir}/${g}_L001_mapq${mapq}.bam
	samtools view -b -q ${mapq} -@${nC} ${bamDir}/${g}_L002_aln.sam > ${bamDir}/${g}_L002_mapq${mapq}.bam
	rm ${bamDir}/${g}_L001_aln.sam ${bamDir}/${g}_L002_aln.sam
done < ${workDir}/gidWells.txt
