#!/usr/bin/bash
#SBATCH -J GATK.HTC --out logs/GATK_HTC_gatk4.%a.log --ntasks 24 --nodes 1 --mem 24G -p intel
#SBATCH --time 1-0:0:0

module unload java
module load java/8
module load gatk/4
module load bcftools/1.9
module load samtools/1.9
module load picard

MEM=32g

ALNFOLDER=aln
HTCFORMAT=cram #default but may switch back to bam
HTCFOLDER=cram # default
HTCEXT=cram
if [ -f config.txt ]; then
    source config.txt
fi
GVCFFOLDER=gvcf_gatk4
DICT=$(echo $REFGENOME | sed 's/fasta$/dict/')

if [ ! -f $DICT ]; then
	picard CreateSequenceDictionary R=$GENOMEIDX O=$DICT
fi
mkdir -p $GVCFFOLDER
TEMP=/scratch
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then 
 echo "need to provide a number by --array slurm or on the cmdline"
 exit
fi

MAX=$(ls $HTCFOLDER/*.$HTCEXT | wc -l | awk '{print $1}')

if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPLESINFO"
 exit
fi
hostname
ALNFILE=$(ls $HTCFOLDER/*.$HTCEXT | sed -n ${N}p)
echo "ALNFILE=$ALNFILE"
if [[ $ALNFILE == "" ]]; then
    echo "cannot find samples in the folder $HTCFOLDER/*.$HTCEXT, exiting ($N)"
    exit
fi
SAMPLE=$(basename $ALNFILE .$HTCEXT)

if [ ! -e $ALNFILE ]; then
    echo "Cannot find $ALNFILE"
    exit
fi
if [ ! -f $GVCFFOLDER/$SAMPLE.g.vcf.gz ]; then
    if [ ! -f $GVCFFOLDER/$SAMPLE.g.vcf ]; then
	# force HTCCALLER to work on BAM files due to some bugs
	# in cram handling by GATK?
	if [ $HTCFORMAT == "bam" ]; then
	    BAM=$ALNFILE
	else
	    BAM=$TEMP/$(basename $ALNFILE .$HTCEXT)".bam"
		time samtools view -O bam --threads $CPU \
	    	--reference $REFGENOME -o $BAM $ALNFILE
	fi    	    
	if [ ! -f $BAM.bai ]; then
		samtools index $BAM
	fi
	time gatk HaplotypeCaller --input $BAM -O $GVCFFOLDER/$SAMPLE.g.vcf \
		--reference $REFGENOME --sample-ploidy 1 \
		-ERC GVCF --sequence-dictionary $DICT \
		 --native-pair-hmm-threads $CPU
	unlink $BAM
    fi
    if [ -f $GVCFFOLDER/$SAMPLE.g.vcf ]; then
	bgzip $GVCFFOLDER/$SAMPLE.g.vcf
	tabix $GVCFFOLDER/$SAMPLE.g.vcf.gz
    fi
fi
