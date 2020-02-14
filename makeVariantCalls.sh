cd /workdir/ns722/
nC=40
mapq=20
# genome="/local/rdata/dbgap_robbins/alfalfa3a.2/3a.2-bionano+unscaffolded-assembly.fa"
# using just the first 978 scaffolds, as they appear to be the main genome sequence. Wrote Andrew Farmer email to confirm
genome="/local/rdata/dbgap_robbins/alfalfaMainGenome/3a.2-bionano+unscaffolded-assembly_1to978.fa"
workDir=alfalfaBioinfo/alfalfaWorkdir
seqDir=alfalfaNovaSeq
bamDir=alfalfaBamsMerged/combinedReps

variantDir=alfalfaVariantCalls

if [ ! -d "$variantDir" ]; then
	mkdir -p $variantDir
fi

ls $bamDir | sed 's|^|'$bamDir'\/|'> ${workDir}/bamNames.txt
sed -i '/\.bai/d' ${workDir}/bamNames.txt

cat ${workDir}/bamNames.txt
wc -l ${workDir}/bamNames.txt

bcftools mpileup --threads ${nC} -a AD -f $genome -b ${workDir}/bamNames.txt  > ${variantDir}/alignAll3a.2.raw.vcf
bcftools call -Ov --threads ${nC} -mv ${variantDir}/alignAll3a.2.raw.vcf > ${variantDir}/alignAll3a.2.call.vcf

