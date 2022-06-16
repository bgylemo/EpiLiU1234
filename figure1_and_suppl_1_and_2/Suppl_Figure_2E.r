library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)

df_reference <- fread("/media/colne37/hippopotamus/thymodevel/data/GTEX/featurecount/res_table.tsv")

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


# source Antonios plot themes
source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")

# read in metadata
metadata <- fread("/media/colne37/hippopotamus/thymodevel/data/GTEX/metadata/METADATA_RNAseq_affymetrix_matched.tsv")

unique(metadata$V1)

#ggplot(metadata, aes(x=V3,y=..count.., fill = V5))+geom_bar(stat="count", position = "dodge") + geom_text(stat="count" ,aes(label=..count..))

################################################################
######## reference-alignment (affymetrix, preprocessed) ########
################################################################
df_affymetrix <- fread("/media/colne37/hippopotamus/thymodevel/data/GTEX/affymetrix_data/GSE45878_MA_data.txt", header = T)

# melt and change colnames
melt_df_affymetrix <- melt(df_affymetrix)
colnames(melt_df_affymetrix) <- c("id_ref", "sample", "val")

# read in affymetrix metadata
meta_affy <- dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/data/GTEX/affymetrix_data/metadata_GSE45878.tsv", header = F), V1, V4, V5, V8)

meta_affy[meta_affy$V8 == "Adipose - Subcutaneous"]$V8 <- "Adipose_Subcutaneous"
meta_affy[meta_affy$V8 == "Artery - Tibial"]$V8 <- "Artery_Tibial"
meta_affy[meta_affy$V8 == "Brain - Cortex"]$V8 <- "Brain_Cortex"
meta_affy[meta_affy$V8 == "Heart - Left Ventricle"]$V8 <- "Heart_Left_Ventricle"
meta_affy[meta_affy$V8 == "Lung"]$V8 <- "Lung"
meta_affy[meta_affy$V8 == "Muscle - Skeletal"]$V8 <- "Muscle_Skeletal"
meta_affy[meta_affy$V8 == "Nerve - Tibial"]$V8 <- "Nerve_Tibial"
meta_affy[meta_affy$V8 == "Skin - Sun Exposed (Lower leg)"]$V8 <- "Skin_Sun_Exposed_Lower_leg"
meta_affy[meta_affy$V8 == "Whole Blood"]$V8 <- "Whole_Blood"

# merge metadata and data, change colnames
melt_df_affymetrix <- merge(melt_df_affymetrix, meta_affy, by.x = c("sample"), by.y = "V1", all = T)
colnames(melt_df_affymetrix) <- c("sample", "id_ref", "val", "sex", "age", "tissue_id")

unique(melt_df_affymetrix[melt_df_affymetrix$tissue_id %in% keep,])

# split
melt_df_affymetrix <- melt_df_affymetrix %>% separate(id_ref, into = c("genes", NA), sep = "_", remove = F)

ensembl_id_to_symbol <- fread("/media/colne37/hippopotamus/thymodevel/data/GTEX/affymetrix_data/ensembl_to_symbol.tsv", header = F)
colnames(ensembl_id_to_symbol) <- c("enseml_id", "gene")
ensembl_id_to_symbol_fix <- ensembl_id_to_symbol %>% separate(enseml_id, into = c("genes", NA), sep = "\\.", remove = F)

melt_df_affymetrix <- merge(melt_df_affymetrix, ensembl_id_to_symbol_fix, by.x =c("genes"), by.y = "genes", all.x =T)

unique(melt_df_affymetrix$gene)


colnames(melt_df_affymetrix) <- c("ensembl_short", "sample", "id_ref", "val", "sex", "age", "tissue_id","ensembl_id", "target_id")


melt_df_affymetrix_short <- dplyr::select(melt_df_affymetrix, target_id, sample, val, sex, tissue_id)


transformation_before_log2FC_calc <- melt_df_affymetrix_short %>% dplyr::group_by(target_id, sex, tissue_id) %>% dplyr::summarise(meanz = mean(val))

sex_split <- dcast(transformation_before_log2FC_calc[transformation_before_log2FC_calc$tissue_id %in% keep,], formula = target_id + tissue_id ~ sex, value.var = "meanz")

sex_split$b <- log2(sex_split$Female / sex_split$Male)

df_MA_sex_short <- dplyr::select(sex_split, -Female, -Male)

# Read in the Tukiainen log2FC data, transform it to long format.
tuki_df <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/Suppl.Table.2.csv")
tuki_df[tuki_df$`Gene name` == "6-Sep"]$`Gene name` <- "SEPT6"
tuki_log2fc_df <- tuki_df %>% dplyr:: select('Gene name', ends_with("logFC"))
melt_tuki_log2fc_df <- melt(tuki_log2fc_df)
melt_tuki_log2fc_df$variable <- gsub(melt_tuki_log2fc_df$variable, pattern = "_logFC", replacement = "")
names(melt_tuki_log2fc_df)[names(melt_tuki_log2fc_df) == 'Gene name'] <- 'gene'

#read in annotations (PAR and the Tukiainen annotation)
PAR <- rbind(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/group-715_PAR1.csv"), fread("/media/colne37/hippopotamus/thymodevel/annotation_data/group-716_PAR2.csv"))
esc <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[,-1]
esc <- esc[!(esc$Gene_name == "IDS" & esc$Reported_XCI_status == "Unknown"),]

# Read in GCF gene annotation
annotation_gene_table <- subset(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv"), type == "gene" & is.na(pseudo))

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"
chrx_annotation[chrx_annotation$Gene_name == "XG",]$category <- "PAR"

chrx_annotation[chrx_annotation$category == "Unknown",]$category <- "Potential"
chrx_annotation[chrx_annotation$category == "Variable",]$category <- "Potential"

################################################################
###################### Add annotation to df ####################
################################################################

df_smash_anno_chrx <- merge(df_MA_sex_short, chrx_annotation, by.x = "target_id", by.y = "Gene_name", all.x = F)

# small fixes
df_smash_anno_chrx$b <- as.numeric(df_smash_anno_chrx$b)

df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Nerve_Tibial",]$tissue_id <- "nerve"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Artery_Tibial",]$tissue_id <- "artery-3"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Heart_Left_Ventricle",]$tissue_id <- "heart-2"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Whole_Blood",]$tissue_id <- "whole blood"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-2"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Adipose_Subcutaneous",]$tissue_id <- "adipose-1"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Brain_Cortex",]$tissue_id <- "brain-6"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Muscle_Skeletal",]$tissue_id <- "muscle"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Lung",]$tissue_id <- "lung"

color_valzz <- c("brain-1" = "#0072B5FF",
                 "brain-2" ="#0072B5FF" ,
                 "brain-3" ="#0072B5FF",
                 "brain-4" = "#0072B5FF", 
                 "brain-5" = "#0072B5FF" ,
                 "brain-6" = "#0072B5FF" ,
                 "brain-7" = "#0072B5FF" ,
                 "brain-8" = "#0072B5FF" ,
                 "brain-8" = "#0072B5FF",
                 "brain-10" = "#0072B5FF",
                 "brain-11" = "#0072B5FF",
                 "brain-12" = "#0072B5FF",
                 "brain-13" ="#0072B5FF",
                 "esophagus-1" = "#4DBBD5FF" ,
                 "esophagus-2" ="#4DBBD5FF",
                 "esophagus-3" = "#4DBBD5FF",
                 "skin-1" ="#00A087FF",
                 "skin-2" = "#00A087FF",
                 "heart-1" = "#3C5488FF",
                 "heart-2" = "#3C5488FF",
                 "adipose-1" = "#F39B7FFF",
                 "adipose-2" ="#F39B7FFF",
                 "artery-1" = "#91D1C2FF",
                 "artery-2" = "#91D1C2FF" ,
                 "artery-3" = "#91D1C2FF",
                 "adrenal gland" = "#8491B4FF",
                 "kidney-1" = "#7E6148FF",
                 "kidney-2" = "#7E6148FF",
                 "bladder" = "#DC0000FF",
                 "nerve" = "#B09C85FF",
                 "salivary gland" = "#BC3C29FF",
                 "throid" = "#E64B35FF",
                 "pituitary" = "#E18727FF",
                 "pancreas" = "#20854EFF",
                 "spleen" = "#7876B1FF",
                 "liver" = "#6F99ADFF",
                 "lung" = "#FFDC91FF",
                 "breast" ="#EE4C97FF",
                 "muscle" = "#374E55FF",
                 "small intestine" = "#DF8F44FF",
                 "colon-1" = "#DF8F44FF",
                 "colon-2" = "#DF8F44FF",
                 "stomach" ="#DF8F44FF",
                 "lymphocytes" = "grey",
                 "fibroblasts" = "#B24745FF",
                 "whole blood" = "#79AF97FF",
                 "thymocytes" = "red")

orderz <- df_smash_anno_chrx %>% dplyr::group_by(tissue_id, category) %>% dplyr::summarise(medianz = median(b))

PAR_orderz <- orderz[orderz$category == "PAR",]
PAR_orderz <- PAR_orderz[order(PAR_orderz$medianz, decreasing = T),]

Escape_orderz <- orderz[orderz$category == "Escape",]
Escape_orderz <- Escape_orderz[order(Escape_orderz$medianz, decreasing = T),]

Inactive_orderz <- orderz[orderz$category == "Inactive",]
Inactive_orderz <- Inactive_orderz[order(Inactive_orderz$medianz, decreasing = T),]

# plot
PAR <- 
  ggplot(df_smash_anno_chrx[df_smash_anno_chrx$category %in% c("PAR"),], aes(x=factor(tissue_id, levels = PAR_orderz$tissue_id), y=b, col = tissue_id)) + 
  geom_quasirandom(show.legend = F, alpha = .5)+
  geom_boxplot(show.legend = F, outlier.shape = NA) + 
  geom_hline(yintercept = 0, lty=2)+
  theme_AL_box_rotX() +
  scale_y_continuous(breaks=c(-0.5,0,0.5), limits = c(-0.51,0.51))+
  scale_color_manual(values = color_valzz)+
  labs(x="", title = "PAR")

escapez <- 
  ggplot(df_smash_anno_chrx[df_smash_anno_chrx$category %in% c("Escape"),], aes(x=factor(tissue_id, levels = Escape_orderz$tissue_id), y=b, col = tissue_id)) + 
  geom_quasirandom(show.legend = F, alpha = .5)+
  geom_boxplot(show.legend = F, outlier.shape = NA) + 
  geom_hline(yintercept = 0, lty=2)+
  theme_AL_box_rotX() +
  scale_y_continuous(breaks=c(-0.5,0,0.5), limits = c(-0.51,0.51))+
  scale_color_manual(values = color_valzz)+
  labs(x="", title = "Escape")

inacz <- 
  ggplot(df_smash_anno_chrx[df_smash_anno_chrx$category %in% c("Inactive"),], aes(x=factor(tissue_id, levels = Inactive_orderz$tissue_id), y=b, col = tissue_id)) + 
  geom_quasirandom(show.legend = F, alpha = .5)+
  geom_boxplot(show.legend = F, outlier.shape = NA) + 
  geom_hline(yintercept = 0, lty=2)+
  theme_AL_box_rotX() +
  scale_y_continuous(breaks=c(-0.5,0,0.5), limits = c(-0.51,0.51))+
  scale_color_manual(values = color_valzz)+
  labs(x="", title = "Inactive")

ggsave2("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_2E.pdf", width = 16, height = 6,
        plot_grid(PAR, escapez, inacz, ncol = 3)
)