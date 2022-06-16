################# Read in and prepare data ####################
# Load libraries
library(data.table)
library(karyoploteR)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplotify)
library(cowplot)

source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")

anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)
anno <- anno[anno$CHR_hg38 == "chrX",]
anno <- dplyr::select(anno,CHR_hg38, Start_hg38, End_hg38, probeID)

anno$Start_hg38 <- as.numeric(anno$Start_hg38)
anno$End_hg38 <- as.numeric(anno$End_hg38)

mefff <- merge(melt(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/femme_only_methylation_beta_values_normalized.tsv")), anno, by = "probeID")


mefff$cell <- "ETP"
mefff[grepl(mefff$variable, pattern = "DN"),]$cell <- "DN"
mefff[grepl(mefff$variable, pattern = "DE"),]$cell <- "DPearly"
mefff[grepl(mefff$variable, pattern = "DL"),]$cell <- "DPlate"
mefff[grepl(mefff$variable, pattern = "S4"),]$cell <- "CD4SP"
mefff[grepl(mefff$variable, pattern = "S8"),]$cell <- "CD8SP"

mefff$sample <- "F1"
mefff[grepl(mefff$variable, pattern = "5"),]$sample <- "F2"
mefff[grepl(mefff$variable, pattern = "6"),]$sample <- "F3"

meff_chrx_meanz <- mefff %>% dplyr::group_by(CHR_hg38, Start_hg38, End_hg38, cell) %>% dplyr::summarise(meanz = mean(value))

ETP_gr <- DataTrack(makeGRangesFromDataFrame(dplyr::select(meff_chrx_meanz[meff_chrx_meanz$cell == "ETP",], -cell), keep.extra.columns=T, ignore.strand=T, seqnames.field = "CHR_hg38", start.field = "Start_hg38", end.field = "End_hg38"), name = "ETP")

DN_gr <- DataTrack(makeGRangesFromDataFrame(dplyr::select(meff_chrx_meanz[meff_chrx_meanz$cell == "DN",], -cell), keep.extra.columns=T, ignore.strand=T, seqnames.field = "CHR_hg38", start.field = "Start_hg38", end.field = "End_hg38"), name = "DN")

DPearly_gr <- DataTrack(makeGRangesFromDataFrame(dplyr::select(meff_chrx_meanz[meff_chrx_meanz$cell == "DPearly",], -cell), keep.extra.columns=T, ignore.strand=T, seqnames.field = "CHR_hg38", start.field = "Start_hg38", end.field = "End_hg38"), name = "DPearly")

DPlate_gr <- DataTrack(makeGRangesFromDataFrame(dplyr::select(meff_chrx_meanz[meff_chrx_meanz$cell == "DPlate",], -cell), keep.extra.columns=T, ignore.strand=T, seqnames.field = "CHR_hg38", start.field = "Start_hg38", end.field = "End_hg38"), name = "DPlate")

CD4SP_gr <- DataTrack(makeGRangesFromDataFrame(dplyr::select(meff_chrx_meanz[meff_chrx_meanz$cell == "CD4SP",], -cell), keep.extra.columns=T, ignore.strand=T, seqnames.field = "CHR_hg38", start.field = "Start_hg38", end.field = "End_hg38"), name = "CD4SP")

CD8SP_gr <- DataTrack(makeGRangesFromDataFrame(dplyr::select(meff_chrx_meanz[meff_chrx_meanz$cell == "CD8SP",], -cell), keep.extra.columns=T, ignore.strand=T, seqnames.field = "CHR_hg38", start.field = "Start_hg38", end.field = "End_hg38"), name = "CD8SP")


itrack <- IdeogramTrack(genome = "hg38", chromosome = "chrX")
gtrack <- GenomeAxisTrack(littleTicks = F)

pdf(file = "/media/colne37/hippopotamus/thymodevel/plots/new_sfig4.pdf", height = 6, width = 6)
plotTracks(list(ETP_gr, DN_gr, DPearly_gr, DPlate_gr, CD4SP_gr, CD8SP_gr, gtrack, itrack), type = "h", chromosome = "chrX", labelPos = "below", exponent = 4)
dev.off()

df_beta <- data.frame(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/femme_only_methylation_beta_values_normalized.tsv", header = T))

# melt
melt_df_beta <- melt(df_beta)

melt_df_beta$cell <- "ETP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DN"),]$cell <- "DN"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DE"),]$cell <- "DPearly"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DL"),]$cell <- "DPlate"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "S4"),]$cell <- "CD4SP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "S8"),]$cell <- "CD8SP"

melt_df_beta$sample <- "F1"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "5"),]$sample <- "F2"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "6"),]$sample <- "F3"

melt_df_beta_chrx <- melt_df_beta[melt_df_beta$probeID %in% anno_short_chrx$probeID,]

melt_df_beta_chrx_meanz <- melt_df_beta_chrx %>% dplyr::group_by(probeID, cell) %>% dplyr::summarise(meanz = mean(value))

mat <- dcast(melt_df_beta_chrx_meanz, probeID ~ cell, value.var = "meanz")
rownames(mat) <- mat$probeID
mat$probeID <- NULL

res <- cor(mat)

library(corrplot)



ggsave2("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_4A_left.pdf", width = 24,
        plot_grid(karyplotz))

pdf("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_4A_right.pdf")
        corrplot.mixed(res, order = 'AOE')
dev.off()