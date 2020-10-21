##This is the main R parsing script for the W7 project 
##It preforms basic data cleaning, visulization, and statistics.
#The main input files are:
#1) The output of snpEff run with TE removal step (W7.snpEff.matrix_NO_TEs.tab).
#2) A literature curated list of putative allergen genes (A_fumigatus_allergen_genes.txt)

#load modules
require(data.table)
require(splitstackshape)
library(plyr)
library(dplyr)
library(tidyverse)
library(rlist)
library(gdata)

#set dir
#setwd()

options(stringsAsFactors = FALSE)
#load main dat files from snpEff
#snpEff<-read.delim("W7.snpEff.matrix.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
snpEff<-read.delim("W7.snpEff.matrix_NO_TEs.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)


###clean the data for easy processing###
#remove the brackets form ALT
snpEff$ALT<-gsub("\\[|\\]", "", snpEff$ALT)

#remove the "/" from all the the SNP calls 
snpEff[] <- lapply(snpEff, gsub, pattern='/', replacement="")

#subset: if the ALT = CEA10, then it's probably not interesting 
#note this will also remove cases where w7 matches the reference (AF293)
#W7 alt not in cea10
to_OG<-snpEff[snpEff$ALT != snpEff$CEA10.2,]
#cea10 alt not in w7
to_OG_invert<-snpEff[snpEff$ALT != snpEff$W72310,]

dim(snpEff)
dim(to_OG)
dim(to_OG_invert)


#get number of no-calls
sum(to_OG$CEA10.2 == ".")
#874
sum(to_OG$W72310 == ".")
#0


#Print the cleaned file for sharing:
write.table(to_OG, "W7_snpEff_cleaned_NO_TEs.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###Investigate the data###
#separate CEA10 from W7 and get totals for mut type 
W7_diff_from_ref<-snpEff[snpEff$REF != snpEff$W72310,]
CEA10_diff_from_ref<-snpEff[snpEff$REF != snpEff$CEA10.2,]

dim(W7_diff_from_ref)
#there are 56036 SNPs in W7 that differ from the reference
dim(CEA10_diff_from_ref)
#there are 13014 SNPs in CEA10 that differ from the reference

#save for table 
n_variants_W7<- nrow(W7_diff_from_ref)
n_variants_CEA10<- nrow(CEA10_diff_from_ref)
summary_table<- data.frame("W72310" = n_variants_W7, "CEA10" = n_variants_CEA10)
row.names(summary_table)<- "Variants"
summary_table<- t(summary_table)


sum(W7_diff_from_ref$W72310 ==".")
#differences from reference include 2454 no-calls in W7
sum(CEA10_diff_from_ref$CEA10.2 ==".")
#differences from reference include 989 no-calls in CEA10

#save for table 
n_nocalls_W7<- sum(W7_diff_from_ref$W72310 ==".")
n_nocalls_CEA10<- sum(CEA10_diff_from_ref$CEA10.2 ==".")
no_calls<- c(n_nocalls_W7,n_nocalls_CEA10)
summary_table<- cbind(summary_table, "no calls" = no_calls)

#call table on each separately 
W7_mut_type<- table(W7_diff_from_ref$TYPE)
W7_mut_type<- data.table(W7_mut_type)
W7_mut_type<- W7_mut_type[order(W7_mut_type$N, decreasing = T),]
dim(W7_mut_type)

CEA10_mut_type<- table(CEA10_diff_from_ref$TYPE)
CEA10_mut_type<- data.table(CEA10_mut_type)
CEA10_mut_type<- CEA10_mut_type[order(CEA10_mut_type$N, decreasing = T),]
dim(CEA10_mut_type)

#call table together
df_w7<- data.frame(cbind(rep("W7", nrow(W7_diff_from_ref)), W7_diff_from_ref$TYPE))
colnames(df_w7)<- "isolate"
df_cea10<- data.frame(cbind(rep("CEA_10", nrow(CEA10_diff_from_ref)), CEA10_diff_from_ref$TYPE))
colnames(df_cea10)<- "isolate"

#combine them
df_all<- rbind(df_w7, df_cea10)
df_all_table<- data.frame(table(df_all))
#df_all_table<- df_all_table[order(df_all_table$Freq, decreasing = T),]

#look at the clades individually 
barplot(df_all_table$Freq,beside=T,ylim=c(0,31000),
        main="mutation type",
        ylab="frequency",
        axis.lty="solid", 
        col = c("#7F5A77", "#E17369"), 
        names.arg = df_all_table$NA.,
        axisnames = TRUE, 
        cex.names = .4, 
        las =2)

legend("topright", 
       legend = c("CEA10", "W7"), 
       fill = c("#7F5A77", "#E17369"))


##Missessence mutations likely have the most impact. 
#look only at missence mutations 
missence_only_W7<- W7_diff_from_ref[W7_diff_from_ref$TYPE == "missense_variant",]
missence_only_CEA10<- CEA10_diff_from_ref[CEA10_diff_from_ref$TYPE == "missense_variant",]


#get mutations in W7 that are not in CEA10
muts_in_W7_and_not_in_CEA10<-W7_diff_from_ref[W7_diff_from_ref$W72310 != W7_diff_from_ref$CEA10.2,]
muts_in_W7_and_not_in_CEA10_clean<- muts_in_W7_and_not_in_CEA10[muts_in_W7_and_not_in_CEA10$W72310 != ".", ]
dim(muts_in_W7_and_not_in_CEA10_clean)

#the inverse:
muts_in_CEA10_and_not_in_W7<-CEA10_diff_from_ref[CEA10_diff_from_ref$CEA10.2 != CEA10_diff_from_ref$W72310,]
muts_in_CEA10_and_not_in_W7_clean<- muts_in_CEA10_and_not_in_W7[muts_in_CEA10_and_not_in_W7$CEA10.2 != ".", ]
dim(muts_in_CEA10_and_not_in_W7_clean)

#save for table 
unique_muts_W7<- nrow(muts_in_W7_and_not_in_CEA10_clean)
unique_muts_CEA10<- nrow(muts_in_CEA10_and_not_in_W7_clean)
unique_calls<- c(unique_muts_W7,unique_muts_CEA10)
summary_table<- cbind(summary_table, "Unique Variants" = unique_calls)

#remove synonymous and intergenic variants 
muts_in_W7_and_not_in_CEA10_clean_a<- muts_in_W7_and_not_in_CEA10_clean[!muts_in_W7_and_not_in_CEA10_clean$TYPE == "intergenic",]
muts_in_W7_and_not_in_CEA10_clean_b<- muts_in_W7_and_not_in_CEA10_clean_a[!muts_in_W7_and_not_in_CEA10_clean_a$TYPE == "intron_variant",]
muts_in_W7_and_not_in_CEA10_no_syn<- muts_in_W7_and_not_in_CEA10_clean_b[!muts_in_W7_and_not_in_CEA10_clean_b$TYPE == "synonymous_variant",]
dim(muts_in_W7_and_not_in_CEA10_no_syn)

#the inverse
muts_in_CEA10_and_not_in_W7_clean_a<- muts_in_CEA10_and_not_in_W7_clean[!muts_in_CEA10_and_not_in_W7_clean$TYPE == "intergenic",]
muts_in_CEA10_and_not_in_W7_clean_b<- muts_in_CEA10_and_not_in_W7_clean_a[!muts_in_CEA10_and_not_in_W7_clean_a$TYPE == "intron_variant",]
muts_in_CEA10_and_not_in_W7_no_syn<- muts_in_CEA10_and_not_in_W7_clean_b[!muts_in_CEA10_and_not_in_W7_clean_b$TYPE == "synonymous_variant",]
dim(muts_in_CEA10_and_not_in_W7_no_syn)

#save for table 
no_syn_muts_W7<- nrow(muts_in_W7_and_not_in_CEA10_no_syn)
no_syn_muts_CEA10<- nrow(muts_in_CEA10_and_not_in_W7_no_syn)
no_syn_calls<- c(no_syn_muts_W7,no_syn_muts_CEA10)
summary_table<- cbind(summary_table, "Unique non-synonamous" = no_syn_calls)

#call table 
muts_in_W7_and_not_in_CEA10_table<- data.frame(table(muts_in_W7_and_not_in_CEA10_no_syn$TYPE))
muts_in_W7_and_not_in_CEA10_table_ordered<-muts_in_W7_and_not_in_CEA10_table[order(muts_in_W7_and_not_in_CEA10_table$Freq, decreasing = T),]
#View(muts_in_W7_and_not_in_CEA10_table_ordered)

#graph these 
barplot(muts_in_W7_and_not_in_CEA10_table_ordered$Freq,beside=T,ylim=c(0,7000),
        main="mutations in W7, not in CEA10 - (no syn. or intergenic mutations)",
        ylab="frequency",
        axis.lty="solid", 
        col = c("#E17369"), 
        names.arg = muts_in_W7_and_not_in_CEA10_table_ordered$Var1,
        axisnames = TRUE, 
        cex.names = .4, 
        las =2)

#load gene list of putative allergen genes
allergen_df<-read.delim("A_fumigatus_allergen_genes.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = F)
allergen_list<- allergen_df$gene_accession

#subset to only those muts that fall into genes in the list of putative allergens
in_allergen_genes<- muts_in_W7_and_not_in_CEA10_no_syn[muts_in_W7_and_not_in_CEA10_no_syn$GENE %in% allergen_df$gene_accession,]
dim(in_allergen_genes)
length(unique(in_allergen_genes$GENE))

in_allergen_genes_CEA10<- muts_in_CEA10_and_not_in_W7_no_syn[muts_in_CEA10_and_not_in_W7_no_syn$GENE %in% allergen_df$gene_accession,]
dim(in_allergen_genes_CEA10)
length(unique(in_allergen_genes_CEA10$GENE))

#print them to txt file 
write.table(in_allergen_genes, "mutations_in_allergen_genes_NO_TEs.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#save for table 
muts_in_allergen_W7<- nrow(in_allergen_genes)
muts_in_allergen_CEA10<- nrow(in_allergen_genes_CEA10)
muts_in_allergen<- c(muts_in_allergen_W7, muts_in_allergen_CEA10)
summary_table<- cbind(summary_table, "Variants in allergen genes" = muts_in_allergen)

#save for table 
in_allergen_W7<- length(unique(in_allergen_genes$GENE))
in_allergen_CEA10<- length(unique(in_allergen_genes_CEA10$GENE))
in_allergen<- c(in_allergen_W7, in_allergen_CEA10)
summary_table<- cbind(summary_table, "Allergen genes with mutations" = in_allergen)


#get the most highly mutated W7 genes 
most_mutations<- data.frame(table(in_allergen_genes$GENE))
most_mutations_ordered<- most_mutations[order(most_mutations$Freq, decreasing = T),]
colnames(most_mutations_ordered)<- c("gene", "mutation_frequency")


#get the most highly mutated CEA10 genes 
most_mutations_CEA10<- data.frame(table(in_allergen_genes_CEA10$GENE))
most_mutations_ordered_CEA10<- most_mutations_CEA10[order(most_mutations_CEA10$Freq, decreasing = T),]
colnames(most_mutations_ordered_CEA10)<- c("gene", "mutation_frequency")

length(intersect(most_mutations_ordered_CEA10$gene, most_mutations_ordered$gene))
intersect(most_mutations_ordered_CEA10$gene, most_mutations_ordered$gene)

# exclusive to W7? yes- the mutations are specific to W7 (not in CEA10 - it doesn't mean that there are NO mutations in this gene in CEA10)
#View(most_mutations_ordered) - add if they these gens are not mutated at all in CEA10. 
#write this for Jane's table
write.table(most_mutations_ordered, "W7_mutated_allergen_genes.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)

#write the summary table 
write.table(summary_table, "W7_summary_table.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)


#get the gene variants in W7 and not in CEA10 most likely to have high impact
high_impact_only<- (muts_in_W7_and_not_in_CEA10_no_syn[muts_in_W7_and_not_in_CEA10_no_syn$IMPACT == "HIGH",])

#how many are high impact?
nrow(high_impact_only)
#in how many unique genes?
length(unique(high_impact_only$GENE))

#write it 
write.table(high_impact_only, "W7_high_impact_only.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)

#write to a file for easy upload inot fungiDB
#high_impact_only_unique_genes<- unique(high_impact_only$GENE)
#write.table(high_impact_only_unique_genes, "W7_high_impact_only_unique_genes.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)
#fungi DB used to generate GO terms for each gene

GO_terms_high_impact<-read.delim("GO_terms_high_impact.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
#View(GO_terms_high_impact)

#search for GO terms: 
#GO: 0006979 (response to oxidative stress)
#GO:0005618 (cell wall)
#GO: GO:0044238 (primary metabolic process)
#melenin 
#GO:0006582 (melanin metabolic process) #(parent to below)
#GO:0042438 (melanin biosynthetic process)


ox_stress_in_high_impact<- GO_terms_high_impact[with(GO_terms_high_impact, grepl("GO:0006979", paste(Computed.GO.Component.IDs, Computed.GO.Function.IDs, Computed.GO.Process.IDs, Curated.GO.Component.IDs, Curated.GO.Function.IDs, Curated.GO.Process.IDs))),]
#View(ox_stress_in_high_impact)
#2

cell_wall_in_high_impact<- GO_terms_high_impact[with(GO_terms_high_impact, grepl("GO:0005618", paste(Computed.GO.Component.IDs, Computed.GO.Function.IDs, Computed.GO.Process.IDs, Curated.GO.Component.IDs, Curated.GO.Function.IDs, Curated.GO.Process.IDs))),]
#View(cell_wall_in_high_impact)
#none

primary_met_in_high_impact<- GO_terms_high_impact[with(GO_terms_high_impact, grepl("GO:0044238", paste(Computed.GO.Component.IDs, Computed.GO.Function.IDs, Computed.GO.Process.IDs, Curated.GO.Component.IDs, Curated.GO.Function.IDs, Curated.GO.Process.IDs))),]
#View(primary_met_in_high_impact)
#none


###do the same for the large list (not just high impact)
all_unique_genes<- unique(muts_in_W7_and_not_in_CEA10_no_syn$GENE)
length(all_unique_genes)
#write.table(all_unique_genes, "W7_all.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)

GO_terms_all<-read.delim("GO_terms_all.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
ox_stress_in_all<- GO_terms_all[with(GO_terms_all, grepl("GO:0006979", paste(Computed.GO.Component.IDs, Computed.GO.Function.IDs, Computed.GO.Process.IDs, Curated.GO.Component.IDs, Curated.GO.Function.IDs, Curated.GO.Process.IDs))),]
dim(ox_stress_in_all)
#20
ox_stress_in_all$Gene.ID

cell_wall_in_all<- GO_terms_all[with(GO_terms_all, grepl("GO:0005618", paste(Computed.GO.Component.IDs, Computed.GO.Function.IDs, Computed.GO.Process.IDs, Curated.GO.Component.IDs, Curated.GO.Function.IDs, Curated.GO.Process.IDs))),]
dim(cell_wall_in_all)
cell_wall_in_all$Gene.ID
#6

primary_met_in_all<- GO_terms_all[with(GO_terms_all, grepl("GO:0044238", paste(Computed.GO.Component.IDs, Computed.GO.Function.IDs, Computed.GO.Process.IDs, Curated.GO.Component.IDs, Curated.GO.Function.IDs, Curated.GO.Process.IDs))),]
dim(primary_met_in_all)
#none

#melanine 
melanin_metabolic_process<- GO_terms_all[with(GO_terms_all, grepl("GO:0006582", paste(Computed.GO.Component.IDs, Computed.GO.Function.IDs, Computed.GO.Process.IDs, Curated.GO.Component.IDs, Curated.GO.Function.IDs, Curated.GO.Process.IDs))),]
dim(melanin_metabolic_process)
#1
melanin_metabolic_process$Gene.ID

melanin_biosynthetic_process<- GO_terms_all[with(GO_terms_all, grepl("GO:0042438", paste(Computed.GO.Component.IDs, Computed.GO.Function.IDs, Computed.GO.Process.IDs, Curated.GO.Component.IDs, Curated.GO.Function.IDs, Curated.GO.Process.IDs))),]
dim(melanin_biosynthetic_process)
#5
melanin_biosynthetic_process$Gene.ID

#print off O2 stress and cell wall genes: 
write.table(ox_stress_in_all, "W7_Ox_stress_GO.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(cell_wall_in_all, "W7_cell_wall_GO.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)

#parse the snpEff file to capture the genes identified above:
ox_stress_in_all_snp<- muts_in_W7_and_not_in_CEA10_no_syn[muts_in_W7_and_not_in_CEA10_no_syn$GENE %in% ox_stress_in_all$Gene.ID,]
cell_wall_in_all_snp<- muts_in_W7_and_not_in_CEA10_no_syn[muts_in_W7_and_not_in_CEA10_no_syn$GENE %in% cell_wall_in_all$Gene.ID,]
mel_in_all_snp<- muts_in_W7_and_not_in_CEA10_no_syn[muts_in_W7_and_not_in_CEA10_no_syn$GENE %in% melanin_biosynthetic_process$Gene.ID,]
mel_in_all_snp2<- muts_in_W7_and_not_in_CEA10_no_syn[muts_in_W7_and_not_in_CEA10_no_syn$GENE %in% melanin_metabolic_process$Gene.ID,]

dim(mel_in_all_snp)
dim(mel_in_all_snp2)

mel_in_all_snp_all<- rbind(mel_in_all_snp, mel_in_all_snp2)

unique(mel_in_all_snp_all$GENE)

#write SNP eff files for these genes:
write.table(ox_stress_in_all_snp, "W7_Ox_stress_snpEff.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(cell_wall_in_all_snp, "W7_cell_wall_wnpEff.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(mel_in_all_snp_all, "W7_mel_snpEff.csv", sep=",", row.names = TRUE, col.names = NA, quote = FALSE)

#get n mutations per gene for Jane's table
ox_stress_in_all_snp_table<- data.frame(table(ox_stress_in_all_snp$GENE))
ox_stress_in_all_snp_table_ordered<- ox_stress_in_all_snp_table[order(ox_stress_in_all_snp_table$Freq, decreasing = T),]
colnames(ox_stress_in_all_snp_table_ordered)<- c("gene", "mutation_frequency")

cell_wall_in_all_snp_table<- data.frame(table(cell_wall_in_all_snp$GENE))
cell_wall_in_all_snp_table_ordered<- cell_wall_in_all_snp_table[order(cell_wall_in_all_snp_table$Freq, decreasing = T),]
colnames(cell_wall_in_all_snp_table_ordered)<- c("gene", "mutation_frequency")

#melanin
#bind metabolic process and bioeynthetic process results
melanin_process<- rbind(melanin_metabolic_process, melanin_biosynthetic_process)
#get totals (n nutations per gene)
mell_in_all_snp_table<- data.frame(table(mel_in_all_snp_all$GENE))
mell_in_all_snp_table_ordered<- mell_in_all_snp_table[order(mell_in_all_snp_table$Freq, decreasing = T),]
colnames(mell_in_all_snp_table_ordered)<- c("gene", "mutation_frequency")

#match up the info with the annotations:
count_info<-read.delim("count_info.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
annotation_info<-read.delim("W7_gene_annotations.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#attach annotation_info$Product.Description
merge_annotations <- left_join(count_info, annotation_info, by = 'gene_ID')

##for mel
#subset to get the correct collums only 
melanin_process_unique<- distinct(melanin_process)
melanin_process_anno<- as.data.frame(cbind(gene = melanin_process_unique$Gene.ID, Product.Description = melanin_process_unique$Product.Description))
#merge annotations with counts
merge_annotations_mel<- left_join(mell_in_all_snp_table_ordered, melanin_process_anno, by = 'gene')

#look for overlap
overlaps<- data.frame(table(merge_annotations$gene_ID))
overlaps_ordered<- overlaps[order(overlaps$Freq, decreasing = T),]

##look at rodA mutations (Afu5g09580)
Afu5g09580<- muts_in_W7_and_not_in_CEA10_no_syn[muts_in_W7_and_not_in_CEA10_no_syn$GENE == "Afu5g09580",]
Afu5g09580_syn<- muts_in_W7_and_not_in_CEA10[muts_in_W7_and_not_in_CEA10$GENE == "Afu5g09580",]

write.table(Afu5g09580, "Afu5g09580_no_syn.csv", sep=",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(Afu5g09580_syn, "Afu5g09580_syn_and_intergenic.csv", sep=",", row.names = FALSE, col.names = TRUE, quote = FALSE)


