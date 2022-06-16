setwd("/media/colne37/seahorse/GTEX/")
library(data.table)
library(dplyr)
library(tidyr)

# source Antonios plot themes
source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")

# read in metadata
metadata <- fread("/media/colne37/hippopotamus/thymodevel/data/GTEX/metadata/metadata_RNAseq_full_for_tissues.tsv")
###################################################################################################################################
############################################################ LOAD DATA ############################################################
###################################################################################################################################

####################################################
######## Pseudoalignment (Salmon -> Sleuth) ########
####################################################
df_pseduo <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/remake/sleuth_remake.tsv")

# remove NAs and headers left while merging
df_pseduo <- df_pseduo[!is.na(df_pseduo$b),]
df_pseduo <- subset(df_pseduo, target_id != "target_id")
df_pseduo$b <- as.numeric(df_pseduo$b)

# make short for merging (i.e. remove unwanted columns)
df_pseduo_short <- dplyr::select(df_pseduo, target_id, b, tissue_id)

#####################################################################
######## reference-alignment (STAR -> featureCount -> edgeR) ########
#####################################################################
df_reference <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/edgeR/edgeR_res_table.tsv")

# make short for merging (i.e. remove unwanted columns) and melt.
df_reference_short <- melt(dplyr::select(df_reference, -logCPM, -'F', -PValue, -enseml_id))


# remove genes with several ensembl_ids, i.e. count above than 1 from this code.
excludez <- df_reference_short %>% dplyr::group_by(gene, variable) %>% dplyr::count()
excludez <- excludez[excludez$n > 1,]

df_reference_short_filt <- df_reference_short[!df_reference_short$gene %in% excludez$gene,]

# change colnames to match pseudo df
colnames(df_reference_short_filt) <- c("target_id", "tissue_id", "b")

# declare vector of tissues to keep
keep <- c(gsub(sort(unique(df_reference_short_filt$tissue_id)), pattern = "logFC.", replacement = ""), "Whole_Blood")


################################################################
####################### Bind together dfs ######################
################################################################

df_pseduo_short$study <- "Pseudo-alignment"
df_reference_short_filt$study <- "Reference-alignment"

df_smash <- rbind(df_pseduo_short, df_reference_short_filt)
df_smash <- df_smash[!is.na(df_smash$target_id),]

df_smash$tissue_id <- gsub(df_smash$tissue_id, pattern = "logFC.", replacement = "")

################################################################
###################### Read in annotations #####################
################################################################

# Read in the Tukiainen log2FC data, transform it to long format.
tuki_df <- fread("/media/colne37/seahorse/GTEX/misc_data/Suppl.Table.2.csv")
tuki_df[tuki_df$`Gene name` == "6-Sep"]$`Gene name` <- "SEPT6"
tuki_log2fc_df <- tuki_df %>% dplyr:: select('Gene name', ends_with("logFC"))
melt_tuki_log2fc_df <- melt(tuki_log2fc_df)
melt_tuki_log2fc_df$variable <- gsub(melt_tuki_log2fc_df$variable, pattern = "_logFC", replacement = "")
names(melt_tuki_log2fc_df)[names(melt_tuki_log2fc_df) == 'Gene name'] <- 'gene'

#read in annotations (PAR and the Tukiainen annotation)
PAR <- rbind(fread("/media/colne37/seahorse/GTEX/indexes_and_annotations/R_annotation/group-715_PAR1.csv"), fread("/media/colne37/seahorse/GTEX/indexes_and_annotations/R_annotation/group-716_PAR2.csv"))
esc <- fread("/media/colne37/seahorse/GTEX/indexes_and_annotations/R_annotation/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[,-1]
esc <- esc[!(esc$Gene_name == "IDS" & esc$Reported_XCI_status == "Unknown"),]

# Read in GCF gene annotation
annotation_gene_table <- subset(fread("/media/colne37/seahorse/GTEX/indexes_and_annotations/R_annotation/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv"), type == "gene" & is.na(pseudo))

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"
chrx_annotation[chrx_annotation$Gene_name == "XG",]$category <- "PAR"

chrx_annotation[chrx_annotation$category == "Unknown",]$category <- "Potential"
#chrx_annotation[chrx_annotation$category == "Variable",]$category <- "Potential"

################################################################
###################### Add annotation to df ####################
################################################################

df_smash_anno_chrx <- merge(df_smash, chrx_annotation, by.x = "target_id", by.y = "Gene_name", all.x = F)
df_smash_anno_all_chrs <- merge(df_smash, chrx_annotation, by.x = "target_id", by.y = "Gene_name", all.x = T)

# small fixes
df_smash_anno_chrx$b <- as.numeric(df_smash_anno_chrx$b)
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Whole_blood",]$tissue_id <- "Whole_Blood"


################################################################
######################### Add TPMs to df #######################
################################################################

# filter metadata

TPMz <- melt(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/remake/tissues_TPM_matrix.tsv"))
metadata_tpm <- dplyr::select(metadata, -2)

metadata_tpm$V2 <- gsub(metadata_tpm$`entity:sample_id`, pattern = "\\.", replacement = "-")

tpmz_mit_meta <- merge(TPMz, metadata_tpm, by.x = "variable", by.y = "entity:sample_id")

tpmz_per_tissue <- tpmz_mit_meta %>% dplyr::group_by(V1, tissue_id) %>% dplyr::summarise(meanz = mean(value))

df_smash_anno_chrx_tpmz <- merge(df_smash_anno_chrx, tpmz_per_tissue, by.x = c("target_id", "tissue_id"), by.y = c("V1", "tissue_id"))

###################################################################################################################################
############################################################ PLOT DATA ############################################################
###################################################################################################################################

library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(ggsci)

# filter metadata
metadata <- metadata[metadata$`entity:sample_id` %in% tpmz_mit_meta$variable,]

#### plot sample counts etc
sample_count_per_sex <- ggplot(metadata, aes(x=sex,y=..count.., fill = factor(sex, levels = c("Male", "Female")))) + geom_bar(stat="count") + geom_text(stat="count", aes(label=..count..)) + theme_AL_box_rotX() + theme(axis.text.x = element_text(size=6), legend.title = element_blank()) + scale_fill_npg()

metadata_tissue_order <- metadata %>% dplyr::group_by(tissue_id) %>% dplyr::count()
metadata_tissue_order <- metadata_tissue_order[order(metadata_tissue_order$n, decreasing = T),]

sample_count_per_sex_and_tissue <- ggplot(metadata, aes(x=factor(tissue_id, levels= metadata_tissue_order$tissue_id),y=..count.., fill = factor(sex, levels = c("Male", "Female")))) + geom_bar(stat="count") + geom_text(stat="count", aes(label=..count..)) + theme_AL_box_rotX() + theme(axis.text.x = element_text(size=6), legend.title = element_blank()) + scale_fill_npg() + labs(x="")


# make stats df
df_smash_anno_chrx_stats <- df_smash_anno_chrx_tpmz[df_smash_anno_chrx_tpmz$meanz > 1,] %>% dplyr::group_by(tissue_id, study, category) %>% rstatix::get_summary_stats(b, type = "common")


# quite a variance in normality, going with median and boxplots
library(rstatix)
options(scipen = 10000)
shappurooo <- df_smash_anno_chrx %>% dplyr::group_by(study, category, tissue_id) %>% shapiro_test(b)





all_cats_are_beautiful <- df_smash_anno_chrx_tpmz[df_smash_anno_chrx_tpmz$meanz > 1,]
all_cats_are_beautiful[all_cats_are_beautiful$b > 1,]$b <- 1.1
all_cats_are_beautiful[all_cats_are_beautiful$b < -1,]$b <- -1.1

allez <- 
  ggplot(all_cats_are_beautiful[all_cats_are_beautiful$study != "Affymetrix",], 
         aes(x=factor(category, levels = c("PAR", "Escape", "Variable", "Inactive", "Potential")),y=b,col=study)) + 
  geom_quasirandom(dodge.width = .75) +  
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept = 0) + 
  theme_AL_box() +
  labs(x="", y="FC (log2)", title = "All tissues")+
  theme(legend.title = element_blank(), legend.position = "top")+
  coord_cartesian(ylim=c(-1.1,1.1))+
  scale_y_continuous(breaks = c(-1,0,1))+
  geom_hline(yintercept = c(-1,1),lty=2)+
  scale_color_manual(values=c("Pseudo-alignment"="black","Reference-alignment"="#5F7D7B"))


ggsave2("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_2B.pdf", height = 8, width = 8,
        plot_grid(allez))

#ggsave2("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_2B_sample_counts.pdf", height = 8, width = 8,        plot_grid(sample_count_per_sex_and_tissue))




