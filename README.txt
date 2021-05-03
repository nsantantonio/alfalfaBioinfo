## minimal workflow

- makeBams.sh (initial sequence alignments to tetraploid genome from Noble with bwa) -> 
- mergeLaneBams.sh (merge bam files across technical and biological replicates with samtools) ->
- makeVariantCalls.sh (call variants with bcftools) - >
- runFilter.sh (initial filter based on read counts, with filterVcf.py) ->
- buildCov.R (apply second filter 0.1 < maf < 0.9 and calaulate Fst, population covariances, R script is messy)
