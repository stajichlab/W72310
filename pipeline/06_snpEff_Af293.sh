#!/usr/bin/bash 
#SBATCH --mem=16G --nodes 1 --ntasks 2 --out logs/snpEff.log -p short
module unload miniconda2
module load miniconda3
module load snpEff
module load bcftools/1.9
module load tabix
source activate pyvcf

SNPEFFOUT=snpEff
SNPEFFGENOME=AfumigatusAf293_FungiDB_39
snpEffConfig=snpEff.config
GFFGENOME=FungiDB-39_AfumigatusAf293.gff
MEM=16g

# this module defines SNPEFFJAR and SNPEFFDIR
if [ -f config.txt ]; then
	source config.txt
fi
GFFGENOMEFILE=$GENOMEFOLDER/$GFFGENOME
FASTAGENOMEFILE=$GENOMEFOLDER/$GENOMEFASTA
if [ -z $SNPEFFJAR ]; then
 echo "need to defined \$SNPEFFJAR in module or config.txt"
 exit
fi
if [ -z $SNPEFFDIR ]; then
 echo "need to defined \$SNPEFFDIR in module or config.txt"
 exit
fi
# could make this a confi

if [ -z $FINALVCF ]; then
	echo "need a FINALVCF variable in config.txt"
	exit
fi

mkdir -p $SNPEFFOUT
if [ ! -e $SNPEFFOUT/$snpEffConfig ]; then
	rsync -a $SNPEFFDIR/snpEff.config $SNPEFFOUT/$snpEffConfig
	echo "# AfumAf293.fungidb " >> $SNPEFFOUT/$snpEffConfig
  	echo "$SNPEFFGENOME.genome : Aspergillus fumigatus Af293 FungiDB" >> $SNPEFFOUT/$snpEffConfig
	chroms=$(grep '##sequence-region' $GFFGENOMEFILE | awk '{print $2}' | perl -p -e 's/\n/, /' | perl -p -e 's/,\s+$/\n/')
	echo -e "\t$SNPEFFGENOME.chromosomes: $chroms" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.mito_A_fumigatus_Af293.codonTable : Mold_Mitochondrial" >> $SNPEFFOUT/$snpEffConfig
	mkdir -p $SNPEFFOUT/data/$SNPEFFGENOME
	gzip -c $GFFGENOMEFILE > $SNPEFFOUT/data/$SNPEFFGENOME/genes.gff.gz
	rsync -aL $REFGENOME $SNPEFFOUT/data/$SNPEFFGENOME/sequences.fa

	java -Xmx$MEM -jar $SNPEFFJAR build -datadir `pwd`/$SNPEFFOUT/data -c $SNPEFFOUT/$snpEffConfig -gff3 -v $SNPEFFGENOME
fi
pushd $SNPEFFOUT
#COMBVCF="../$FINALVCF/$PREFIX.selected_nofixed.SNP.vcf.gz ../$FINALVCF/$PREFIX.selected_nofixed.INDEL.vcf.gz"
COMBVCF="../$FINALVCF/$PREFIX.selected.SNP.vcf.gz ../$FINALVCF/$PREFIX.selected.INDEL.vcf.gz"
for n in $COMBVCF
do
 echo $n
 st=$(echo $n | perl -p -e 's/\.gz//')
 if [ ! -f $n ]; then
	 bgzip $st
	 tabix $n
 fi
done
INVCF=$PREFIX.comb_selected.SNP.vcf
OUTVCF=$PREFIX.snpEff.vcf
OUTTAB=$PREFIX.snpEff.tab
OUTDOMAIN=$PREFIX.snpEff_protdomain.tsv
bcftools concat -a -d both -o $INVCF -O v $COMBVCF

java -Xmx$MEM -jar $SNPEFFJAR eff -dataDir `pwd`/data -v $SNPEFFGENOME $INVCF > $OUTVCF

bcftools view -Oz -o $OUTVCF.gz $OUTVCF
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT]\t%INFO/ANN\n' $OUTVCF.gz > $OUTTAB

echo "running matrix dump"
../../scripts/snpEff_2_tab.py $OUTVCF ../$REFGENOME > $PREFIX.snpEff.matrix.tab
pigz -f -k $PREFIX.snpEff.matrix.tab $OUTTAB
../../scripts/map_snpEff2domains.py --vcf $OUTVCF --domains ../genome/FungiDB-39_AfumigatusAf293_InterproDomains.txt --output $OUTDOMAIN
