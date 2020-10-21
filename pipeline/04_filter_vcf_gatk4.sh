#!/usr/bin/bash
#SBATCH --nodes 1 -p short
#SBATCH --ntasks 1
#SBATCH --mem 16G
#SBATCH --job-name=GATK.select_filter
#SBATCH --output=logs/select_filter_gatk4.log

module load gatk/4
module unload java
module load java/8
module load bcftools/1.9

hostname

CONFIG=config.txt
if [ -f $CONFIG ]; then
    source $CONFIG
fi
GENOME=$GENOMEFOLDER/$GENOMEFASTA
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
    CPU=1
fi

if [ -z $FINALVCF ]; then
	echo "Need FINALVCF for output"
	exit
fi
echo "$PREFIX $INFILE"
PREFIX=$PREFIX
INFILE=$FINALVCF/$PREFIX.all.vcf.gz
INSNP=$FINALVCF/$PREFIX.SNP.vcf
ININDEL=$FINALVCF/$PREFIX.INDEL.vcf
FILTEREDSNP=$FINALVCF/$PREFIX.filtered.SNP.vcf
FILTEREDINDEL=$FINALVCF/$PREFIX.filtered.INDEL.vcf
SNPONLY=$FINALVCF/$PREFIX.selected.SNP.vcf
INDELONLY=$FINALVCF/$PREFIX.selected.INDEL.vcf

FILTEREDFIXEDSNP=$FINALVCF/$PREFIX.filtered_fixed.SNP.vcf
FILTEREDFIXEDINDEL=$FINALVCF/$PREFIX.filtered_fixed.INDEL.vcf
SNPNOFIXED=$FINALVCF/$PREFIX.selected_nofixed.SNP.vcf
INDELNOFIXED=$FINALVCF/$PREFIX.selected_nofixed.INDEL.vcf


 if [ ! -f $INSNP.gz ]; then
  gatk SelectVariants \
  -R $REFGENOME \
  --variant $INFILE \
  -O $INSNP \
  --restrict-alleles-to BIALLELIC \
  --select-type-to-include SNP
  bgzip $INSNP
  tabix $INSNP.gz
 fi
 if [ ! -f $ININDEL.gz ]; then
  gatk SelectVariants \
  -R $REFGENOME \
  --variant $INFILE \
  --output $ININDEL \
  --select-type-to-include INDEL --select-type-to-include MIXED --select-type-to-include MNP
  bgzip $ININDEL
  tabix $ININDEL.gz
 fi

 if [ ! -f $FILTEREDSNP.gz ]; then
   gatk VariantFiltration --output $FILTEREDSNP \
   --variant $INSNP.gz -R $REFGENOME \
   --cluster-window-size 10  \
   --filter-expression "QD < 2.0" --filter-name QualByDepth \
   --filter-expression "MQ < 40.0" --filter-name MapQual \
   --filter-expression "QUAL < 100" --filter-name QScore \
   --filter-expression "MQRankSum < -12.5" --filter-name MapQualityRankSum \
   --filter-expression "SOR > 3.0" --filter-name StrandOddsRatio \
   --filter-expression "FS > 60.0" --filter-name FisherStrandBias \
   --filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRank \
   --missing-values-evaluate-as-failing
   bgzip $FILTEREDSNP
   tabix $FILTEREDSNP.gz
 fi

 if [ ! -f $FILTEREDFIXEDSNP.gz ]; then
   gatk VariantFiltration --output $FILTEREDFIXEDSNP \
   --variant $INSNP.gz -R $REFGENOME \
   --cluster-window-size 10  \
   --filter-expression "AF > 0.99" --filter-name FixedAllele \
   --filter-expression "QD < 2.0" --filter-name QualByDepth \
   --filter-expression "MQ < 40.0" --filter-name MapQual \
   --filter-expression "QUAL < 100" --filter-name QScore \
   --filter-expression "MQRankSum < -12.5" --filter-name MapQualityRankSum \
   --filter-expression "SOR > 3.0" --filter-name StrandOddsRatio \
   --filter-expression "FS > 60.0" --filter-name FisherStrandBias \
   --filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRank \
   --missing-values-evaluate-as-failing
   bgzip $FILTEREDFIXEDSNP
   tabix $FILTEREDFIXEDSNP.gz

 fi
 if [ ! -f $FILTEREDINDEL.gz ]; then
  gatk VariantFiltration --output $FILTEREDINDEL \
  --variant $ININDEL.gz -R $REFGENOME \
  --cluster-window-size 10  -filter "QD < 2.0" --filter-name QualByDepth \
  -filter "MQRankSum < -12.5" --filter-name MapQualityRankSum \
  -filter "SOR > 4.0" --filter-name StrandOddsRatio \
  -filter "FS > 200.0" --filter-name FisherStrandBias \
  -filter "ReadPosRankSum < -20.0" --filter-name ReadPosRank 
  # -filter "InbreedingCoeff < -0.8" --filter-name InbreedCoef \

  bgzip $FILTEREDINDEL
  tabix $FILTEREDINDEL.gz
 fi

 if [ ! -f $FILTEREDFIXEDINDEL.gz ]; then
  gatk VariantFiltration --output $FILTEREDFIXEDINDEL \
  --variant $ININDEL.gz -R $REFGENOME \
   --filter-expression "AF > 0.99" --filter-name FixedAllele \
  --cluster-window-size 10  -filter "QD < 2.0" --filter-name QualByDepth \
  -filter "MQRankSum < -12.5" --filter-name MapQualityRankSum \
  -filter "SOR > 4.0" --filter-name StrandOddsRatio \
  -filter "FS > 200.0" --filter-name FisherStrandBias \
  -filter "ReadPosRankSum < -20.0" --filter-name ReadPosRank 

  #-filter "InbreedingCoeff < -0.8" --filter-name InbreedCoef \
  bgzip $FILTEREDFIXEDINDEL
  tabix $FILTEREDFIXEDINDEL.gz
 fi

 if [ ! -f $SNPONLY.gz ]; then
   gatk SelectVariants -R $REFGENOME \
   --variant $FILTEREDSNP.gz \
   --output $SNPONLY \
   --exclude-filtered
   bgzip $SNPONLY
   tabix $SNPONLY.gz
 fi

 if [ ! -f $INDELONLY.gz ]; then
   gatk SelectVariants -R $REFGENOME \
   --variant $FILTEREDINDEL.gz \
   --output $INDELONLY \
   --exclude-filtered 
   bgzip $INDELONLY
   tabix $INDELONLY.gz
 fi

 if [ ! -f $SNPNOFIXED.gz ]; then
     gatk SelectVariants -R $REFGENOME \
	 --variant $FILTEREDFIXEDSNP.gz \
	 --output $SNPNOFIXED \
	 --exclude-filtered
    bgzip $SNPNOFIXED
    tabix $SNPNOFIXED.gz
 fi

 if [ ! -f $INDELNOFIXED.gz ]; then
     gatk SelectVariants -R $REFGENOME \
	 --variant $FILTEREDFIXEDINDEL.gz \
	 --output $INDELNOFIXED \
	 --exclude-filtered
     bgzip $INDELNOFIXED
     tabix $INDELNOFIXED.gz
 fi
