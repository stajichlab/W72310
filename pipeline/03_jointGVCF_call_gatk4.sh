#!/usr/bin/bash
#SBATCH --mem 32G --nodes 1 --ntasks 2 -J GATK.GVCFGeno --out logs/GVCFGenoGATK4_test.log --time 3-0:00:00

MEM=128g
module unload java
module load java/8
module load picard
module load gatk/4
module load tabix
module load parallel

FINALVCF=vcf
mkdir -p $FINALVCF
if [ -f config.txt ]; then
	source config.txt
fi
OUT=$FINALVCF/$PREFIX.all.vcf

if [ ! -f $REFGENOME ]; then
	module load samtools/1.9
	samtools faidx $REFGENOME
fi
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=2
fi
if [[ $(ls $GVCFFOLDER | grep -c -P "\.g.vcf$") -gt "0" ]]; then
	parallel -j $CPU bgzip {} ::: $GVCFFOLDER/*.g.vcf
	parallel -j $CPU tabix -f {} ::: $GVCFFOLDER/*.g.vcf.gz
fi

FILES=$(ls $GVCFFOLDER/*.g.vcf.gz | sort | perl -p -e 's/(\S+)\n/-V $1 /')
INTERVALS=$(cut -f1 $REFGENOME.fai  | perl -p -e 's/(\S+)\n/--intervals $1 /g')
DB=${GVCFFOLDER}_db
if [ ! -f $OUT.gz ]; then
 if [ ! -f $OUT ]; then
	 rm -rf $DB
	 gatk GenomicsDBImport --genomicsdb-workspace-path $DB $FILES $INTERVALS
	 time gatk GenotypeGVCFs --reference $REFGENOME --output $OUT -V gendb://$DB
  else
  	echo "not rerunning GenotypeGVCFs as $OUT already exists"
  fi
  bgzip $OUT
  tabix $OUT.gz
else
	echo "not re-running since $OUT.gz already exists"
fi
