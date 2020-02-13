
nC=40
# bwa index ref.fa

# bwa mem ref.fa reads.fq > aln-se.sam

# bwa mem ref.fa read1.fq read2.fq > aln-pe.sam

# bwa aln ref.fa short_read.fq > aln_sa.sai

# bwa samse ref.fa aln_sa.sai short_read.fq > aln-se.sam

# bwa sampe ref.fa aln_sa1.sai aln_sa2.sai read1.fq read2.fq > aln-pe.sam

# bwa bwasw ref.fa long_read.fq > aln.sam 

# with 20 cpus, this takes 4 hours for one. yikes

# Chilean seed increase 
bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.sam
bwa mem -t 50 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_aln.sam

bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_aln.sam
bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_aln.sam

samtools view -b -q 20 -@40 alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.sam > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.bam
samtools view -b -q 20 -@40 alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_aln.sam > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_mapq20.bam

samtools view -b -q 20 -@40 alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_aln.sam > alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_mapq20.bam
samtools view -b -q 20 -@40 alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_aln.sam > alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_mapq20.bam


samtools merge -@40 alfalfaAlignCADL/Chilean.bam alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L002_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-A02-Chilean-Rep1_S9_L002_mapq20.bam

samtools sort -T /tmp/aln.sorted -@40 -o alfalfaAlignCADL/Chilean_sorted.bam alfalfaAlignCADL/Chilean.bam
bedtools genomecov -ibam alfalfaAlignCADL/Chilean_sorted.bam -bg > alfalfaAlignCADL/Chilean_intervalCoverage.tab # must be sorted first for use of bam

bedtools genomecov -ibam alfalfaAlignCADL/Chilean_sorted.bam -d > alfalfaAlignCADL/Chilean_bpCoverage.tab # 

################

# Chilean original
nC=40

for i in alfalfaDNAext-D01-ChileanOriginal-Rep1_S4 alfalfaDNAext-D07-ChileanOriginal-Rep3_S28; do 
	bwa mem -t ${nC} alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/${i}_L001_R1_001.fastq.gz alfalfaNovaSeq/${i}_L001_R2_001.fastq.gz > alfalfaAlignCADL/${i}_L001_aln.sam
	bwa mem -t ${nC} alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/${i}_L002_R1_001.fastq.gz alfalfaNovaSeq/${i}_L002_R2_001.fastq.gz > alfalfaAlignCADL/${i}_L002_aln.sam

	samtools view -b -q 20 -@${nC} alfalfaAlignCADL/${i}_L001_aln.sam > alfalfaAlignCADL/${i}_L001_mapq20.bam
	samtools view -b -q 20 -@${nC} alfalfaAlignCADL/${i}_L002_aln.sam > alfalfaAlignCADL/${i}_L002_mapq20.bam
done

samtools merge -@${nC} alfalfaAlignCADL/ChileanOriginal.bam alfalfaAlignCADL/alfalfaDNAext-D01-ChileanOriginal-Rep1_S4_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-D01-ChileanOriginal-Rep1_S4_L002_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-D07-ChileanOriginal-Rep3_S28_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-D07-ChileanOriginal-Rep3_S28_L002_mapq20.bam

samtools sort -T /tmp/aln.sorted -@${nC} -o alfalfaAlignCADL/ChileanOriginal_sorted.bam alfalfaAlignCADL/ChileanOriginal.bam



#####

# M falcata original

bwa mem -t ${nC} alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_aln.sam
bwa mem -t ${nC} alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_aln.sam

samtools view -b -q 20 -@${nC} alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_aln.sam > alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_mapq20.bam
samtools view -b -q 20 -@${nC} alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_aln.sam > alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_mapq20.bam

samtools merge -@${nC} alfalfaAlignCADL/MfalcataOriginal.bam alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-C01-MfalcataOriginal-Rep1_S3_L002_mapq20.bam

samtools sort -T /tmp/aln.sorted -@${nC} -o alfalfaAlignCADL/MfalcataOriginal_sorted.bam alfalfaAlignCADL/MfalcataOriginal.bam
bedtools genomecov -ibam alfalfaAlignCADL/MfalcataOriginal_sorted.bam -bg > alfalfaAlignCADL/MfalcataOriginal_intervalCoverage.tab # must be sorted first for use of bam

bedtools genomecov -ibam alfalfaAlignCADL/MfalcataOriginal_sorted.bam -d > alfalfaAlignCADL/MfalcataOriginal_bpCoverage.tab # must be sorted first for use of bam

################


# M falcata seed increase 
nC=20
bwa mem -t ${nC} alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L001_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L001_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L001_aln.sam
bwa mem -t ${nC} alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L002_R1_001.fastq.gz alfalfaNovaSeq/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L002_R2_001.fastq.gz > alfalfaAlignCADL/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L002_aln.sam

samtools view -b -q 20 -@${nC} alfalfaAlignCADL/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L001_aln.sam > alfalfaAlignCADL/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L001_mapq20.bam
samtools view -b -q 20 -@${nC} alfalfaAlignCADL/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L002_aln.sam > alfalfaAlignCADL/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L002_mapq20.bam

samtools merge -@${nC} alfalfaAlignCADL/Mfalcata.bam alfalfaAlignCADL/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-E02-Mfalcata-Rep1_S13_L002_mapq20.bam

samtools sort -T /tmp/aln.sorted -@${nC} -o alfalfaAlignCADL/Mfalcata_sorted.bam alfalfaAlignCADL/Mfalcata.bam
bedtools genomecov -ibam alfalfaAlignCADL/Mfalcata_sorted.bam -bg > alfalfaAlignCADL/Mfalcata_intervalCoverage.tab # must be sorted first for use of bam

bedtools genomecov -ibam alfalfaAlignCADL/Mfalcata_sorted.bam -d > alfalfaAlignCADL/Mfalcata_bpCoverage.tab # must be sorted first for use of bam

################


# ch x mf

nC=40

# alfalfaDNAext-E01-ChileanxMfalcata-Rep1_S5_L001_R1_001.fastq.gz
# alfalfaDNAext-E07-ChileanxMfalcata-Rep3_S29_L001_R1_001.fastq.gz
for i in alfalfaDNAext-E01-ChileanxMfalcata-Rep1_S5 alfalfaDNAext-E07-ChileanxMfalcata-Rep3_S29; do 
	bwa mem -t ${nC} alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/${i}_L001_R1_001.fastq.gz alfalfaNovaSeq/${i}_L001_R2_001.fastq.gz > alfalfaAlignCADL/${i}_L001_aln.sam
	bwa mem -t ${nC} alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaNovaSeq/${i}_L002_R1_001.fastq.gz alfalfaNovaSeq/${i}_L002_R2_001.fastq.gz > alfalfaAlignCADL/${i}_L002_aln.sam

	samtools view -b -q 20 -@${nC} alfalfaAlignCADL/${i}_L001_aln.sam > alfalfaAlignCADL/${i}_L001_mapq20.bam
	samtools view -b -q 20 -@${nC} alfalfaAlignCADL/${i}_L002_aln.sam > alfalfaAlignCADL/${i}_L002_mapq20.bam
done

samtools merge -@${nC} alfalfaAlignCADL/ChileanxMfalcata.bam alfalfaAlignCADL/alfalfaDNAext-E01-ChileanxMfalcata-Rep1_S5_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-E01-ChileanxMfalcata-Rep1_S5_L002_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-E07-ChileanxMfalcata-Rep3_S29_L001_mapq20.bam alfalfaAlignCADL/alfalfaDNAext-E07-ChileanxMfalcata-Rep3_S29_L002_mapq20.bam

samtools sort -T /tmp/aln.sorted -@${nC} -o alfalfaAlignCADL/ChileanxMfalcata_sorted.bam alfalfaAlignCADL/ChileanxMfalcata.bam


####################################
# samtools mpileup -uf alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta  alfalfaAlignCADL/Chilean_sorted.bam alfalfaAlignCADL/MfalcataOriginal_sorted.bam | bcftools call -mv > alfalfaAlignCADL/chmf.raw.vcf
# bcftools filter -s LowQual -e '%QUAL<20 || DP>100' var.raw.vcf  > chmf.flt.vcf
    
# cd alfalfaAlignCADL
samtools index -@30 alfalfaAlignCADL/Chilean_sorted.bam
samtools index -@30 alfalfaAlignCADL/ChileanOriginal_sorted.bam
samtools index -@30 alfalfaAlignCADL/MfalcataOriginal_sorted.bam
samtools index -@30 alfalfaAlignCADL/Mfalcata_sorted.bam
samtools index -@30 alfalfaAlignCADL/ChileanxMfalcata_sorted.bam
# cd ..

bcftools mpileup --threads 50 -a AD -f alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta  alfalfaAlignCADL/Chilean_sorted.bam alfalfaAlignCADL/ChileanOriginal_sorted.bam alfalfaAlignCADL/MfalcataOriginal_sorted.bam alfalfaAlignCADL/Mfalcata_sorted.bam alfalfaAlignCADL/ChileanxMfalcata_sorted.bam > alfalfaAlignCADL/chmf.raw.vcf
bcftools call -Ov --threads 50 -mv > alfalfaAlignCADL/chmf.ploidy.vcf


#####


samtools view -b -h alfalfaAlignCADL/MfalcataOriginal_sorted.bam 'scaffold45|size1037999|quiver' > alfalfaAlignCADL/MFscf45.bam
samtools view -b -h alfalfaAlignCADL/Chilean_sorted.bam 'scaffold45|size1037999|quiver' > alfalfaAlignCADL/CHscf45.bam

bcftools mpileup --threads 30 -a AD -f alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta  alfalfaAlignCADL/CHscf45.bam alfalfaAlignCADL/MFscf45.bam > chmf.scf45.vcf
bcftools call -Ov --threads 30 --ploidy 4 -mv > alfalfaAlignCADL/chmf.scf45.vcf


bcftools mpileup --threads 30 -uf alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta  alfalfaAlignCADL/Chilean_sorted.bam alfalfaAlignCADL/MfalcataOriginal_sorted.bam | bcftools call --threads 30 -mv > alfalfaAlignCADL/chmf.raw.vcf



    # samtools mpileup -uf alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta alfalfaAlignCADL/MfalcataOriginal.bam alfalfaAlignCADL/Chilean.bam | bcftools view -bvcg - > var.raw.bcf  
    # bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf  


################################

bwa mem -t 8 genome.fa reads.fastq | samtools sort -@8 -o output.bam -

samtools view -b -q 20 alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.sam > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.bam
samtools sort -T /tmp/aln.sorted -@20 -o alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20_sorted.bam alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.bam 
bedtools genomecov -ibam alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.bam -bg > alfalfaDNAext-A08-Chilean-Rep3_S33_L001_intervalCoverage.tab # must be sorted first for use of bam
# bedtools genomecov alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq20.sam -bg > alfalfaDNAext-A08-Chilean-Rep3_S33_L001_intervalCoverage.tab

awk '{ print $5 }' alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.sam > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_mapq.txt

# awk -v x=5 '{print $x}' > alfalfaAlignCADL/alfalfaDNAext-A08-Chilean-Rep3_S33_L001_aln.txt 


samtools merge finalBamFile.bam *.bam


 samtools mpileup [-EB] [-C capQcoef] [-r reg] [-f in.fa] [-l list] [-Q minBaseQ] [-q minMapQ] in.bam [in2.bam [...]] 


 bedtools coverage [OPTIONS] -a <FILE> \











for i in 
 # bwa mem -t 40 alfalfaCADLv0.95/CADL_HM342.v0.95P.fasta ${i}.fasta.gz | samtools sort -@40 -o output.bam -
