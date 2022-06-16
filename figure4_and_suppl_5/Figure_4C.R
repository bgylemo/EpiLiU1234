options(stringsAsFactors = F)

# Load data
source("/media/colne37/hippopotamus/thymodevel/my_scripts/figure4_and_suppl_5/combined_Turner_Thymocyte_load_sleuth_FC_NM_AND_NR_sans_chrY_PAR_transcripts.R")

# read in results from normal thymocyte development analysis (see /media/colne37/hippopotamus/thymodevel/my_scripts/figure1_and_suppl_1_and_2/Suppl_Figure_1D_HM_and_Cluster.R)
res_lrt <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/res_lrt.tsv")[,-1]

# Filter out genes not expressed in any subtype
tpm_filter <- txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo$Celltype,function(y) mean(y>=1) )==1) ),]


# average per sex (i.e. ETP F1, F2, F3 into F_ETP and so on)
Femme <- c("Rep1", "Rep4", "Rep5")
Homme <- c("Rep0", "Rep2", "Rep3")

df_tpm_filter <- data.frame(tpm_filter)

turner_columns <- dplyr::select(df_tpm_filter, Turner_ETP, Turner_DN, Turner_DPE, Turner_DPL, Turner_SP4, Turner_SP8)
colnames(turner_columns) <- c("ETP_Turner", "DN_Turner", "DPearly_Turner", "DPlate_Turner", "CD4SP_Turner", "CD8SP_Turner")

ETP_F <- dplyr::select(df_tpm_filter, TD2018_ETP_Rep1, TD2018_ETP_Rep4, TD2018_ETP_Rep5)
DN_F <- dplyr::select(df_tpm_filter, TD2018_DN_Rep1, TD2018_DN_Rep4, TD2018_DN_Rep5)
DPearly_F <- dplyr::select(df_tpm_filter, TD2018_DPearly_Rep1, TD2018_DPearly_Rep4, TD2018_DPearly_Rep5)
DPlate_F <- dplyr::select(df_tpm_filter, TD2018_DPlate_Rep1, TD2018_DPlate_Rep4, TD2018_DPlate_Rep5)
CD4SP_F <- dplyr::select(df_tpm_filter, TD2018_CD4SP_Rep1, TD2018_CD4SP_Rep4, TD2018_CD4SP_Rep5)
CD8SP_F <- dplyr::select(df_tpm_filter, TD2018_CD8SP_Rep1, TD2018_CD8SP_Rep4, TD2018_CD8SP_Rep5)

ETP_F$meanz <- rowMeans(ETP_F)
DN_F$meanz <- rowMeans(DN_F)
DPearly_F$meanz <- rowMeans(DPearly_F)
DPlate_F$meanz <- rowMeans(DPlate_F)
CD4SP_F$meanz <- rowMeans(CD4SP_F)
CD8SP_F$meanz <- rowMeans(CD8SP_F)

femme_means <- dplyr::select(cbind(ETP_F, DN_F$meanz, DPearly_F$meanz, DPlate_F$meanz, CD4SP_F$meanz, CD8SP_F$meanz), -TD2018_ETP_Rep1, -TD2018_ETP_Rep4, -TD2018_ETP_Rep5)
#kek$gene <- rownames(kek)
colnames(femme_means) <- c("ETP_F", "DN_F", "DPearly_F", "DPlate_F", "CD4SP_F", "CD8SP_F")


ETP_M <- dplyr::select(df_tpm_filter, TD2018_ETP_Rep2, TD2018_ETP_Rep3)
DN_M <- dplyr::select(df_tpm_filter, TD2018_DN_Rep2, TD2018_DN_Rep3)
DPearly_M <- dplyr::select(df_tpm_filter, TD2018_DPearly_Rep0, TD2018_DPearly_Rep2, TD2018_DPearly_Rep3)
DPlate_M <- dplyr::select(df_tpm_filter, TD2018_DPlate_Rep0, TD2018_DPlate_Rep2, TD2018_DPlate_Rep3)
CD4SP_M <- dplyr::select(df_tpm_filter, TD2018_CD4SP_Rep0, TD2018_CD4SP_Rep2, TD2018_CD4SP_Rep3)
CD8SP_M <- dplyr::select(df_tpm_filter, TD2018_CD8SP_Rep0, TD2018_CD8SP_Rep2, TD2018_CD8SP_Rep3)

ETP_M$meanz <- rowMeans(ETP_M)
DN_M$meanz <- rowMeans(DN_M)
DPearly_M$meanz <- rowMeans(DPearly_M)
DPlate_M$meanz <- rowMeans(DPlate_M)
CD4SP_M$meanz <- rowMeans(CD4SP_M)
CD8SP_M$meanz <- rowMeans(CD8SP_M)

homme_means <- dplyr::select(cbind(ETP_M, DN_M$meanz, DPearly_M$meanz, DPlate_M$meanz, CD4SP_M$meanz, CD8SP_M$meanz), -TD2018_ETP_Rep2, -TD2018_ETP_Rep3)
#kek$gene <- rownames(kek)
colnames(homme_means) <- c("ETP_M", "DN_M", "DPearly_M", "DPlate_M", "CD4SP_M", "CD8SP_M")

TPM_sex_means <- cbind(femme_means, homme_means, turner_columns)

# Calculate Z-scores
df_z <- t(apply(TPM_sex_means[row.names(TPM_sex_means) %in% res_lrt$target_id,],1,function(x) (x-mean(x))/sd(x) ))

order_colz <- c("ETP_F","ETP_M","ETP_Turner",
                "DN_F","DN_M","DN_Turner",
                "DPearly_F","DPearly_M","DPearly_Turner",
                "DPlate_F","DPlate_M","DPlate_Turner",
                "CD4SP_F","CD4SP_M","CD4SP_Turner",
                "CD8SP_F","CD8SP_M","CD8SP_Turner")

cl <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/clusterz.txt")

test_cl <- cl[cl$V1 %in% rownames(df_z),]

df_z_cut <- df_z[rownames(df_z) %in% cl$V1,]

df_z_cut_ordered <- df_z_cut[order(match(rownames(df_z_cut), test_cl$V1)), , drop = FALSE]


test_cl$x <- factor(test_cl$x, levels=c("1","2","3","4","5","6","7"))

library(ComplexHeatmap)

test_lol <- data.frame(namez = rownames(df_z_cut_ordered), index = 1:nrow(df_z_cut_ordered))

test_lol_kek <- test_lol[test_lol$namez %in% c( "DNMT3B", "HES1", "HOXA9","LYL1", "MYB",  "MYCN",  "NOTCH1",  "SPI1",  "TLX2",  "BCL11A", "DNMT1",  "LMO2",  "RUNX1",  "STAT4",  "BACH2",  "CCR4",  "CCR5",  "CCR7",  "CD28",  "DNMT3A",  "ETS1",  "FOXP3",  "IKZF1",  "IL2RA",  "IRF4","MAF",  "NFATC1",  "NFATC2",  "RUNX3",  "STAT1",  "STAT2",  "STAT3",  "STAT5A",  "STAT5B",  "STAT6",  "TBX21",  "TNF",  "BCL2L1",  "NFATC3", "RAG1",  "RORC",  "TCF12",  "BCL11B",  "CD3E",  "CD3G",  "E4F1",  "GATA3",  "KDM5A",  "LEF1",  "TET1"),]

ha = rowAnnotation(foo = anno_mark(at = test_lol_kek$index, labels = test_lol_kek$namez))

library(circlize)
library(RColorBrewer)

pdf("/media/colne37/hippopotamus/thymodevel/plots/Figure_4C.pdf",width = 12,height = 16)
set.seed(13)
Heatmap(df_z_cut_ordered, split=test_cl$x, show_row_names = F,column_order = order_colz, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))), show_row_dend = F, right_annotation = ha, use_raster = T, cluster_row_slices = FALSE, border = T)
dev.off()
