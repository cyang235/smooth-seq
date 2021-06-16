##################################################################
# For Smooth-seq data analysis
# Scripts are provided for transparency
# Need be modified to fit other datasets
# Updated on Feb 22, 2021
# 
# Copyright (c) 2021, Cheng Yang, TangLab@PKU
# 
# Permission is hereby granted, free of charge, to 
# any person obtaining a copy of this software and 
# associated documentation files (the "Software"), 
# to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software 
# is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice 
# shall be included in all copies or substantial 
# portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH 
# THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
###################################################################


#1. Generate CCS from subreads
# ccs v4.0.0 was used
# https://github.com/PacificBiosciences/ccs 
pbindex ${subread}
ccs $subread ${prefix}.ccs.bam --minPasses 1 -j 0


#2. Demultiplex and trim barcoded CCS reads for each cell
# ${primer} need be provied, according to your samples
lima --dump-clips --same --peek-guess --split-bam-named -j 0 ${input}.ccs.bam $primer ${out}.bam


#3. Align CCS reads to reference genome for each cell
pbindex ${cell_bam}
bam2fastq -o ${cell_fq} ${cell_bam}
pbmm2 align ${ref_genome_index} ${cell_fq}.fastq.gz ${cell_id}.bam --preset CCS --sort -j 20 


#4. SV detection and benchmark
# pbsv v2.2.2 was used 
# https://github.com/PacificBiosciences/pbsv
# pbsv discover stage was run with GRCh38 tandem repeat annotations
# https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
pbsv discover --tandem-repeats ${repeats} ${cell_id}.bam ${cell_id}.svsig.gz
pbsv call --ccs ${reference} ${cell_id}.svsig.gz ${cell_id}.var.vcf


# SV can be filtered with, for example bcftools
# For Smooth-seq, SV calls flagged with IMPRECISE or SHADOWED were filtered
# The length of SV need be longer than 100 bp
# Only PASS calls were considered
# SV calls need to be supported by at least two reads
# e.g.,
bcftools filter -i 'FORMAT/AD[0:1] > 1 & IMPRECISE=0 & SHADOWED=0' ${cell_vcf} > ${cell_filtered_vcf}
# Please use your filter expressions


# merge SVs with SURVIVOR
# e.g.,
SURVIVOR merge fileList 1000 2 1 -1 -1 -1 ${prefix}.vcf


# benchmark SVs with Truvari
# https://github.com/spiralgenetics/truvari
truvari bench -b ${base}.vcf.gz \
    -c  ${cell}.vcf.gz \
    -r 1000 \
    -f ${REF} \
    --multimatch --sizemax 10000000 \
    -o ${cell_id}


#5. SNV and indel detection and benchmark
# GATK HaplotypeCaller v4.1.4.0 was used
# PCRINDELMODEL="AGGRESSIVE"
# SNP="gatk4/bundle/hg38/dbsnp_146.hg38.vcf.gz"
# knownsites1="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
# knownsites2="dbsnp_146.hg38.vcf.gz"
# knownsites3="1000G_phase1.snps.high_confidence.hg38.vcf.gz"


gatk --java-options "-Xmx8g -Xms8g"  HaplotypeCaller \
     --reference "${reference}" \
     --input "${cell_bam}" \
     --output ${cell_id}.g.vcf.gz \
     --sample-ploidy "${PLOIDY}" \
     --pcr-indel-model "${PCRINDELMODEL}" \
     --read-filter MappingQualityReadFilter \
     --read-filter NotSecondaryAlignmentReadFilter \
     --minimum-mapping-quality "${MAPQ}" \
     --emit-ref-confidence GVCF \
     --annotation QualByDepth \
     --annotation RMSMappingQuality \
     --annotation MappingQualityRankSumTest \
     --annotation FisherStrand \
     --annotation StrandOddsRatio \
     --annotation Coverage \
     --annotation ReadPosRankSumTest \
     --annotation-group StandardAnnotation \
     --annotation-group StandardHCAnnotation


gatk --java-options "-Xmx12g -Xms8g" GenotypeGVCFs \
     --output ${cell_id}.vcf.gz \
     --reference ${reference} \
     --dbsnp ${SNP} \
     --sample-ploidy "${PLOIDY}" \
     --standard-min-confidence-threshold-for-calling 2.0 \
     --annotation QualByDepth \
     --annotation RMSMappingQuality \
     --annotation MappingQualityRankSumTest \
     --annotation FisherStrand \
     --annotation StrandOddsRatio \
     --annotation Coverage \
     --annotation ReadPosRankSumTest \
     --variant ${gvcf}


# select SNPs
gatk --java-options "-Xmx12g -Xms8g"  SelectVariants \
     --reference "${REF}" \
     --variant "${cellvcf}" \
     --select-type-to-include SNP \
     --output "${cell_id}.SNP.raw.vcf"

# select 1bp indels
gatk --java-options "-Xmx12g -Xms8g" SelectVariants \
     --reference "${REF}" \
     --variant "${cellvcf}" \
     --select-type-to-include INDEL \
     --max-indel-size 1 \
     --output "${cell_id}.1bpINDEL.raw.vcf"

# select >1bp indels
gatk --java-options "-Xmx12g -Xms8g" SelectVariants \
     --reference "${REF}" \
     --variant "${cellvcf}" \
     --select-type-to-include INDEL \
     --min-indel-size 2 \
     --output "${cell_id}.longINDEL.raw.vcf"

# filter SNPs
gatk --java-options "-Xmx12g -Xms8g"  VariantFiltration \
     --reference "${REF}" \
     --variant "${cell_id}.SNP.raw.vcf" \
     --output "${cell_id}.SNP.filtered.vcf" \
     --filter-name "QDlt2" \
     --filter-expression "QD < 2.0"

# filter 1bp indels
gatk --java-options "-Xmx12g -Xms8g"  VariantFiltration \
     --reference "${REF}" \
     --variant "${cell_id}.1bpINDEL.raw.vcf" \
     --output "${cell_id}.1bpINDEL.filtered.vcf" \
     --filter-name "QDlt5" \
     --filter-expression "QD < 5.0"

# filter >1bp indels
gatk --java-options "-Xmx12g -Xms8g"  VariantFiltration \
     --reference "${REF}" \
     --variant "${cell_id}.longINDEL.raw.vcf" \
     --output "${cell_id}.longINDEL.filtered.vcf" \
     --filter-name "QDlt2" \
     --filter-expression "QD < 2.0"


# benchmark SNV with vcfeval
# https://github.com/RealTimeGenomics/rtg-tools
rtg vcfeval -b ${bulkvcf} \
     -c ${cellvcf}  \
     -o ${cell_id} \
     -t ${REF} \
     --sample ALT,ALT \
     --squash-ploidy
