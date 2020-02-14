cd /workdir/ns722/
nC=40
mapq=20
# genome="/local/rdata/dbgap_robbins/alfalfa3a.2/3a.2-bionano+unscaffolded-assembly.fa"
# using just the first 978 scaffolds, as they appear to be the main genome sequence. Wrote Andrew Farmer email to confirm
genome="/local/rdata/dbgap_robbins/alfalfaMainGenome/3a.2-bionano+unscaffolded-assembly_1to978.fa"
workDir=alfalfaBioinfo/alfalfaWorkdir
seqDir=alfalfaNovaSeq
techRepDir=alfalfaBamsMerged/techReps
combinedDir=alfalfaBamsMerged/combinedReps
bamDir=alfalfaBams
# techRepDir=alfalfaBamsMerged/testtechReps
# combinedDir=alfalfaBamsMerged/testcombinedReps
# bamDir=testBams

sed -i '/BLANK/d' ${workDir}/gidWells.txt
sed -i '/BLANK/d' ${workDir}/gid.txt


if [ ! -d "$techRepDir" ]; then
	mkdir -p $techRepDir
fi

if [ ! -d "$combinedDir" ]; then
	mkdir -p $combinedDir
fi

# g=alfalfaDNAext-A01-AfricanxMfalcata-Rep1_S1
# cat ${workDir}/gidWells.txt

while IFS= read -r g; do
	samtools merge -@${nC} ${techRepDir}/${g}.bam ${bamDir}/${g}_L001_mapq20.bam ${bamDir}/${g}_L002_mapq20.bam 
	samtools sort -T /tmp/${g}aln.sorted -@${nC} -o ${techRepDir}/${g}_sorted.bam ${techRepDir}/${g}.bam
	rm ${techRepDir}/${g}.bam
done < ${workDir}/gidWells.txt


while IFS= read -r h; do
	ls $techRepDir | grep "\\-${h}-Rep1\|\\-${h}-Rep3" | sed 's|^|'$techRepDir'\/|' > tmpNames
	samtools merge -@${nC} ${combinedDir}/${h}.bam -b tmpNames
	samtools sort -T /tmp/${h}aln.sorted -@${nC} -o ${combinedDir}/${h}_sorted.bam ${combinedDir}/${h}.bam
	rm tmpNames
	rm ${combinedDir}/${h}.bam
done < ${workDir}/gid.txt


ls ${combinedDir}/*.bam | xargs -n1 -P${nC} samtools index


# for i in ${combinedDir}/*.bam; do
# 	samtools index -@${nC} ${i}
# done

# ls -lh $techRepDir | head
# ls -lh $techRepDir
# ls $combinedDir


# rm $techRepDir/*
# rm $combinedDir/*

# 	rm ${techRepDir}/${g}.bam
# 	if [ g == "Mfalcata" ]; then
# 		samtools merge -@${nC} ${techRepDir}/${g}.bam ${bamDir}/${g}_L001_mapq20.bam ${bamDir}/${g}_L002_mapq20.bam 
# 	else 
# 		samtools merge -@${nC} ${techRepDir}/${g}.bam ${bamDir}/${g}_L001_mapq20.bam ${bamDir}/${g}_L002_mapq20.bam 
# 	fi

# 	fi
#   echo "gid: $g; wellLane: $gw"

# done <  ${workDir}/gid.txt 3<  ${workDir}/gidWells.txt

# awk '!a[$0]++' ${workDir}/gid.txt > ${workDir}/uniqueGids.txt


# awk '!a[$0]++' ${workDir}/gid.txt > ${workDir}/uniqueGids.txt

# # awk '{!seen[$0]++};END{for(i in seen) if(seen[i]==1)print i}' ${workDir}/gid.txt

# grep ${workDir}/gid.txt

# g=alfalfaDNAext-E02-Mfalcata-Rep1_S13

# sed "s/_.*\|alfalfaDNAext-...-//g" ${workDir}/gidWells.txt

# echo "-${g}-Rep1" 
# grep "\\-Chilean-" ${workDir}/gidWells.txt

# Chilean


# while IFS= read -r g && IFS= read -r gw <&3; do
#   echo "gid: $g; wellLane: $gw"

# done <  ${workDir}/gid.txt 3<  ${workDir}/gidWells.txt