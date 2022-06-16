#Load sleuth object for differential expression
library(data.table)
library(sleuth)
library(dplyr)
library(tidyr)
library(cowplot)

source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")

######## read in data #############
message("loading sleuth object")
so_thymo_turner <- sleuth_load("/media/colne37/hippopotamus/thymodevel/data/Turner/combined_genes_Turner_Thymo_sans_chrY_PAR_NM_AND_NR_ONLY_log2_sleuth")

df_tpm <- melt(sleuth_to_matrix(so_thymo_turner, 'obs_norm', 'tpm'))

colnames(df_tpm) <- c("target_id", "cell_sample", "TPM")

######## read in annotation #############
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

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"
chrx_annotation[chrx_annotation$Gene_name == "XG",]$category <- "PAR"
chrx_annotation[chrx_annotation$category == "Unknown",]$category <- "Potential"
chrx_annotation[chrx_annotation$category == "Variable",]$category <- "Potential"


# keep only chrX genes
df_tpm_filt <- df_tpm[df_tpm$target_id %in% chrx_annotation$Gene_name | df_tpm$target_id == "PUDP",]

# fixes
df_tpm_filt$cell_sample <- gsub(df_tpm_filt$cell_sample, pattern = "TD2018_", replacement = "")
df_tpm_filt <- df_tpm_filt %>% separate(cell_sample, into = c("cell", "sample"), sep = "_")
df_tpm_filt[df_tpm_filt$cell == "Turner" & df_tpm_filt$sample == "ETP",]$cell <- "ETP"
df_tpm_filt[df_tpm_filt$cell == "Turner" & df_tpm_filt$sample == "DN",]$cell <- "DN"
df_tpm_filt[df_tpm_filt$cell == "Turner" & df_tpm_filt$sample == "DPE",]$cell <- "DPearly"
df_tpm_filt[df_tpm_filt$cell == "Turner" & df_tpm_filt$sample == "DPL",]$cell <- "DPlate"
df_tpm_filt[df_tpm_filt$cell == "Turner" & df_tpm_filt$sample == "SP4",]$cell <- "CD4SP"
df_tpm_filt[df_tpm_filt$cell == "Turner" & df_tpm_filt$sample == "SP8",]$cell <- "CD8SP"
df_tpm_filt[df_tpm_filt$sample %in% c("DN","DPE","DPL","ETP","SP4","SP8"),]$sample <- "Turner"

# add sex
df_tpm_filt$sex <- ifelse(df_tpm_filt$sample %in% c("Rep1", "Rep4", "Rep5"), yes = "female", no = "male")
df_tpm_filt[df_tpm_filt$sample == "Turner",]$sex <- "Turner"

df_tpm_filt_anno <- merge(df_tpm_filt, chrx_annotation, by.x = "target_id", by.y = "Gene_name", all.x = T)

# find genes not expressed in any subtype
removez <- df_tpm_filt_anno %>% dplyr::group_by(target_id, cell) %>% dplyr::summarise(meanz = mean(TPM))
removez <- removez[removez$meanz < 1,] %>% dplyr::count() %>% subset(n == 6)

df_tpm_filt_anno <- df_tpm_filt_anno[!df_tpm_filt_anno$target_id %in% removez$target_id,]

df_tpm_filt_anno_fix_pudp_sept6 <- df_tpm_filt_anno
df_tpm_filt_anno_fix_pudp_sept6[df_tpm_filt_anno_fix_pudp_sept6$target_id == "PUDP",]$category <- "Escape"
df_tpm_filt_anno_fix_pudp_sept6[df_tpm_filt_anno_fix_pudp_sept6$target_id == "SEPT6",]$category <- "Escape"

df_tpm_filt_anno_stats <- df_tpm_filt_anno_fix_pudp_sept6 %>% dplyr::group_by(category, cell, sex) %>% rstatix::get_summary_stats(TPM, type = "common")

## PLOT ##
top <- ggplot(df_tpm_filt_anno_stats[df_tpm_filt_anno_stats$category != "Potential" ,], aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP")), y = mean,col = sex)) + geom_pointrange(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(width = .75)) + facet_wrap(~factor(category, levels = c("PAR", "Escape", "Inactive")), ncol = 3) + theme_AL_box_rotX() + labs(x="", y="mean TPM +-SE") + coord_cartesian(ylim=c(0,300))


ggsave2("/media/colne37/hippopotamus/thymodevel/plots/Figure_4A.pdf",height = 4, width = 8,
        plot_grid(top)
)

