#!/bin/bash

# sorting and indexing

samtools sort Control.bam > Control.sorted.bam
samtools index Control.sorted.bam

samtools sort Tumor.bam > Tumor.sorted.bam
samtools index Tumor.sorted.bam

# count reads
samtools view -c Control.sorted.bam
samtools view -c Tumor.sorted.bam

# count reads with quality >30
samtools view -c -q 30 Control.sorted.bam
samtools view -c -q 30 Tumor.sorted.bam

# statistics on reads (avarage quality, insert size std dev, ...)
samtools stats Control.sorted.bam > stat.control.txt; less stat.control.txt
samtools stats Tumor.sorted.bam > stat.tumor.txt; less stat.tumor.txt

# explore coverage 
samtools bedcov Captured_Regions.bed Control.sorted.bam > BEDCov.Control.CR.txt | less -S
samtools bedcov Captured_Regions.bed Tumor.sorted.bam > BEDCov.Tumor.CR.txt | less -S

# realignment

java -jar ../tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.bam -o realigner.intervals.Control -L Captured_Regions.bed
java -jar ../tools/GenomeAnalysisTK.jar -T IndelRealigner -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.bam -targetIntervals Control.realigner.intervals -o Control.sorted.realigned.bam -L Captured_Regions.bed

java -jar ../tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.bam -o realigner.intervals.Tumor -L Captured_Regions.bed
java -jar ../tools/GenomeAnalysisTK.jar -T IndelRealigner -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.bam -targetIntervals Tumor.realigner.intervals -o Tumor.sorted.realigned.bam -L Captured_Regions.bed

# count realigned reads
samtools view Control.sorted.realigned.bam | grep OC | wc -l
samtools view Tumor.sorted.realigned.bam | grep OC | wc -l

# recalibration after realignment

java -jar ../tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I ../Realignment/Control.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -o recal.table.Control -L Captured_Regions.bed
java -jar ../tools/GenomeAnalysisTK.jar -T PrintReads -R ../Annotations/human_g1k_v37.fasta -I ../Realignment/Control.sorted.realigned.bam -BQSR recal.table.Control -o Control.sorted.realigned.recalibrated.bam -L Captured_Regions.bed --emit_original_quals
java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I ../Realignment/Control.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -BQSR recal.table -o after_recal.table.Control -L Captured_Regions.bed
java -jar ../Tools/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ../Annotations/human_g1k_v37.fasta -before recal.table.Control -after after_recal.table.Control -csv recal.Control.csv -plots recal.Control.pdf

java -jar ../tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I ../Realignment/Tumor.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -o recal.table.Tumor -L Captured_Regions.bed
java -jar ../tools/GenomeAnalysisTK.jar -T PrintReads -R ../Annotations/human_g1k_v37.fasta -I ../Realignment/Tumor.sorted.realigned.bam -BQSR Tumor.recal.table -o Tumor.sorted.realigned.recalibrated.bam -L Captured_Regions.bed --emit_original_quals
java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -BQSR Tumor.recal.table -o Tumor.after_recal.table -L Captured_Regions.bed
java -jar ../Tools/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ../Annotations/human_g1k_v37.fasta -before recal.table.Tumor -after after_recal.table.Tumor -csv recal.Tumor.csv -plots recal.Tumor.pdf

# deduplication (picard)

java -jar ../../Tools/picard.jar MarkDuplicates I=Control.sorted.bam O=Control.sorted.dedup.bam REMOVE_DUPLICATES=true TMP_DIR=/tmp METRICS_FILE=Control.picard.log ASSUME_SORTED=true
samtools index Control.sorted.dedup.bam

java -jar ../../Tools/picard.jar MarkDuplicates I=Tumor.sorted.bam O=Tumor.sorted.dedup.bam REMOVE_DUPLICATES=true TMP_DIR=/tmp METRICS_FILE=Tumor.picard.log ASSUME_SORTED=true
samtools index Tumor.sorted.dedup.bam

# somatic copy number calling

samtools mpileup -q 1 -f ../Annotations/human_g1k_v37.fasta Control.sorted.realigned.recalibrated.debup.bam Tumor.sorted.realigned.recalibrated.debup.bam | java - jar ../Tools/VarScan.v2.3.9.jar copynumber --output-file SCNA --mpileup 1

java -jar ../Tools/VarScan.v2.3.9.jar copyCaller SCNA.copynumber --output-file SCNA.copynumber.called 

Rscript CBS.R

# variant calling

bcftools mpileup -Ou -a DP -f ../../Annotations/human_g1k_v37.fasta Control.sorted.bam | bcftools call -Ov -c -v > Control.BCF.vcf 

bcftools mpileup -Ou -a DP -f ../../Annotations/human_g1k_v37.fasta Tumor.sorted.bam | bcftools call -Ov - c -v > Tumor.BCF.vcf 

java -jar ../../Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../../Annotations/human_g1k_v37.fasta -I Control.sorted.bam -o Control.GATK.vcf -L chr20.bed 

java -jar ../../Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../../Annotations/human_g1k_v37.fasta -I Tumor.sorted.bam -o Tumor.GATK.vcf -L chr20.bed 

vcftools --minQ 20 --min-meanDP 30 --remove-indels --vcf Control.BCF.vcf --out Control.BCF --recode --recode INFO-all

vcftools --minQ 20 --min-meanDP 30 --remove-indels --vcf Tumor.BCF.vcf --out Tumor.BCF --recode --recode INFO-all 

vcftools --minQ 20 --min-meanDP 30 --remove-indels --vcf Control.GATK.vcf --out Control.GATK --recode -- recode-INFO-all 

vcftools --minQ 20 --min-meanDP 30 --remove-indels --vcf Tumor.GATK.vcf --out Tumor.GATK --recode -- recode-INFO-all

vcftools --vcf Control.BCF.recode.vcf --diff Control.GATK.recode.vcf --diff-site 
vcftools --vcf Tumor.BCF.recode.vcf --diff Tumor.GATK.recode.vcf --diff-site 
vcftools --vcf Control.BCF.recode.vcf --diff Tumor.BCF.recode.vcf --diff-site
vcftools --vcf Control.GATK.recode.vcf --diff Tumor.GATK.recode.vcf --diff-site

# variant annotation 

java -Xmx4g -jar ../Tools/snpEff/snpEff.jar -v hg19kg ../05_VariantCalling/Data/Control.BCF.recode.vcf - s Control.BCF.recode.ann.html > Control.BCF.recode.ann.vcf 

java -Xmx4g -jar ../Tools/snpEff/snpEff.jar -v hg19kg ../05_VariantCalling/Data/Control.GATK.recode.vcf -s Control.GATK.recode.ann.html > Control.GATK.recode.ann.vcf 

java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar Annotate ../Annotations/hapmap_3.3.b37.vcf Control.BCF.recode.ann.vcf > Control.BCF.recode.ann2.vcf

java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar Annotate ../Annotations/hapmap_3.3.b37.vcf Control.GATK.recode.ann.vcf > Control.GATK.recode.ann2.vcf 

java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar Annotate ../Annotations/clinvar_Pathogenic.vcf Control.BCF.recode.ann2.vcf > Control.BCF.recode.ann3.vcf 

java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar Annotate ../Annotations/clinvar_Pathogenic.vcf Control.GATK.recode.ann2.vcf > Control.GATK.recode.ann3.vcf 

cat Sample.BCF.recode.ann3.vcf | java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar filter "(ANN[ANY].IMPACT = 'HIGH') & (DP > 20) & (exists ID)" 
cat Sample.GATK.recode.ann3.vcf | java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar filter "(ANN[ANY].IMPACT = 'HIGH') & (DP > 20) & (exists ID)" 

cat Sample.BCF.recode.ann3.vcf | java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar filter "(exists CLNSIG)"
cat Sample.GATK.recode.ann3.vcf | java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar filter "(exists CLNSIG)"

# somatic variant calling

samtools mpileup -B -f ../Annota<ons/human_g1k_v37.fasta Control.sorted.realigned.recalibrated.debup.bam > Control.sorted.pileup

java -jar ../Tools/VarScan.v2.3.9.jar mpileup2snp Control.pileup --p-value 0.01 --output-vcf 1 > Control.VARSCAN.vcf

vcftools --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf Control.VARSCAN.vcf --out Control.VARSCAN --recode --recode-INFO-all

vcftools --diff Control.BCF.recode.vcf --vcf Control.VARSCAN.recode.vcf --diff-site 

samtools mpileup -q 1 -f ../Annotations/human_g1k_v37.fasta Tumor.sorted.realigned.recalibrated.debup.bam > Tumor.sorted.pileup 

java -jar ../Tools/VarScan.v2.3.9.jar somatic Control.sorted.pileup Tumor.sorted.pileup --output-snp somatic.pm --output-indel somatic.indel --output-vcf 1 

java -Xmx4g -jar ../../Tools/snpEff/snpEff.jar -v hg19kg somatic.pm.vcf -s somatic.pm.vcf.html > somatic.pm.ann.vcf

# ancestry analysis

Rscript Run.R

# purity ploidy estimation 

bcftools mpileup -Ou -a DP -f ../Annotations/human_g1k_v37.fasta Control.sorted.bam | bcftools call -Ov -c -v > Control.BCF.vcf

bcftools mpileup -Ou -a DP -f ../Annotations/human_g1k_v37.fasta Tumor.sorted.bam | bcftools call -Ov -c -v > Tumor.BCF.vcf 

grep -E "(^#|0/1)" Control.BCF.vcf > Control.het.vcf 
grep -E "(^#|0/1)" Tumor.BCF.vcf > Tumor.het.vcf 

java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter -R ../Annota<ons/human_g1k_v37.fasta -o recal.Control.csv -I Control.sorted.realigned.recalibrated.debup.bam -sites Control.het.GATK.vcf -U ALLOW_N_CIGAR_READS -minDepth 20 --minMappingQuality 20 --minBaseQuality 20 

java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter -R ../Annota<ons/human_g1k_v37.fasta -o Tumor.csv -I Tumor.sorted.realigned.recalibrated.debup.bam -sites Control.het.vcf -U ALLOW_N_CIGAR_READS -minDepth 20 --minMappingQuality 20 --minBaseQuality 20 

Rscript Clonet.R

