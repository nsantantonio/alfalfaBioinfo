# bwa index ref.fa

# bwa mem ref.fa reads.fq > aln-se.sam

# bwa mem ref.fa read1.fq read2.fq > aln-pe.sam

# bwa aln ref.fa short_read.fq > aln_sa.sai

# bwa samse ref.fa aln_sa.sai short_read.fq > aln-se.sam

# bwa sampe ref.fa aln_sa1.sai aln_sa2.sai read1.fq read2.fq > aln-pe.sam

# bwa bwasw ref.fa long_read.fq > aln.sam 

# with 20 cpus, this takes 4 hours for one. yikes
bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.sam
bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_aln.sam

bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_aln.sam
bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_aln.sam

bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_aln.sam
bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_aln.sam

bwa mem -t 8 genome.fa reads.fastq | samtools sort -@8 -o output.bam -

samtools view -b -q 20 alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.sam > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.bam
samtools sort -T /tmp/aln.sorted -o alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20_sorted.bam alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.bam 
bedtools genomecov -ibam alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.bam -bg > alfalfaDNAext-A08-Chilean-Rep3_S33_L001_intervalCoverage.tab # must be sorted first for use of bam
# bedtools genomecov alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.sam -bg > alfalfaDNAext-A08-Chilean-Rep3_S33_L001_intervalCoverage.tab

awk '{ print $5 }' alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.sam > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq.txt

# awk -v x=5 '{print $x}' > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.txt 


samtools merge finalBamFile.bam *.bam


 samtools mpileup [-EB] [-C capQcoef] [-r reg] [-f in.fa] [-l list] [-Q minBaseQ] [-q minMapQ] in.bam [in2.bam [...]] 


 bedtools coverage [OPTIONS] -a <FILE> \











for i in 
 # bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta ${i}.fasta.gz | samtools sort -@40 -o output.bam -
