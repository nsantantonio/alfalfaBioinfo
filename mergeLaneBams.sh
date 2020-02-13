nC=40
mapq=20
# genome="/local/rdata/dbgap_robbins/alfalfa3a.2/3a.2-bionano+unscaffolded-assembly.fa"
# using just the first 978 scaffolds, as they appear to be the main genome sequence. Wrote Andrew Farmer email to confirm
genome="/local/rdata/dbgap_robbins/alfalfaMainGenome/3a.2-bionano+unscaffolded-assembly_1to978.fa"
workDir=alfalfaWorkdir
seqDir=alfalfaNovaSeq
bamDir=alfalfaBamsMerged

while read g; do
	# bwa mem -t ${nC} ${genome} ${seqDir}/${g}_L001_R1_001.fastq.gz ${seqDir}/${g}_L001_R2_001.fastq.gz > ${bamDir}/${g}_L001_aln.sam
	# bwa mem -t ${nC} ${genome} ${seqDir}/${g}_L002_R1_001.fastq.gz ${seqDir}/${g}_L002_R2_001.fastq.gz > ${bamDir}/${g}_L002_aln.sam
	# samtools view -b -q ${mapq} -@${nC} ${bamDir}/${g}_L001_aln.sam > ${bamDir}/${g}_L001_mapq${mapq}.bam
	# samtools view -b -q ${mapq} -@${nC} ${bamDir}/${g}_L002_aln.sam > ${bamDir}/${g}_L002_mapq${mapq}.bam
	# rm ${bamDir}/${g}_L001_aln.sam ${bamDir}/${g}_L002_aln.sam

	samtools merge -@${nC} alfalfaAlignCADL/ChileanOriginal.bam alfalfaAlignCADL/alfalfaDNAext-D01-ChileanOriginal-Rep1_S4_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-D01-ChileanOriginal-Rep1_S4_L002_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-D07-ChileanOriginal-Rep3_S28_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-D07-ChileanOriginal-Rep3_S28_L002_mapq20.bam
	samtools sort -T /tmp/aln.sorted -@${nC} -o alfalfaAlignCADL/ChileanOriginal_sorted.bam alfalfaAlignCADL/ChileanOriginal.bam



done < ${workDir}/gidWells.txt

