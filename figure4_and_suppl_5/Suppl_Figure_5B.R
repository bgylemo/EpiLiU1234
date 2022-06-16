options(stringsAsFactors = F)

source("/media/colne37/hippopotamus/thymodevel/my_scripts/data_prep_figure1/3_Thymocyte_load_sleuth_FC_NM_AND_NR_sans_chrY_PAR_transcripts.R")


sleuth_contrast <- function(so,x,ref){
  new_s2c <- so$sample_to_covariates
  new_s2c$Celltype <- as.factor(new_s2c$Celltype)
  new_s2c$Celltype <- relevel(new_s2c$Celltype, ref = ref)
  
  so2 <- so
  so2$sample_to_covariates <- new_s2c
  new_fit_name <- paste0('ref',ref)
  so2 <- sleuth_fit(so2,~Celltype, fit_name = new_fit_name)
  so2 <- sleuth_wt(so2, paste0('Celltype',x), which_model = new_fit_name)
  return(so2)
}

so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"DN","ETP")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"DPearly","DN")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"DPlate","DPearly")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"CD4SP","DPlate")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"CD8SP","DPlate")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"CD4SP","CD8SP")


res_cntr <- list(
  "DN_vs_ETP" = sleuth_results(so_thymo_transcripts,'CelltypeDN', which_model = "refETP",test_type = "wt",show_all = F, pval_aggregate = F),
  "DPearly_vs_DN" = sleuth_results(so_thymo_transcripts,'CelltypeDPearly', which_model = "refDN",test_type = "wt",show_all = F, pval_aggregate = F),
  "DPlate_vs_DPearly" =sleuth_results(so_thymo_transcripts,'CelltypeDPlate', which_model = "refDPearly",test_type = "wt",show_all = F, pval_aggregate = F),
  "CD4SP_vs_DPlate" = sleuth_results(so_thymo_transcripts,'CelltypeCD4SP', which_model = "refDPlate",test_type = "wt",show_all = F, pval_aggregate = F),
  "CD8SP_vs_DPlate" = sleuth_results(so_thymo_transcripts,'CelltypeCD8SP', which_model = "refDPlate",test_type = "wt",show_all = F, pval_aggregate = F),
  "CD4SP_vs_CD8SP" = sleuth_results(so_thymo_transcripts,'CelltypeCD4SP', which_model = "refCD8SP",test_type = "wt",show_all = F, pval_aggregate = F)
)

names(res_cntr)

##############################
## Differential gene counts ##
##############################

tpm_filter <- txi_thymo_transition$abundance[apply(txi_thymo_transition$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]

tpm_avg <- t(apply(tpm_filter,1,function(x) tapply(x,meta_thymo_transition$Celltype, mean) ))

log2r <- function(res_nm,cutoff=1e-3){
  idx <- unlist(strsplit(res_nm,"_"))
  res <- na.omit(res_cntr[[res_nm]])
  nm <- intersect(res$gene_id[res$qval < cutoff],row.names(tpm_avg))
  x <- tpm_avg[nm,idx[1]]
  y <- tpm_avg[nm,idx[3]]
  # is the gene expressed in any condition?
  isexpr <- x>=1|y>=1
  l2r <- log2(x+1) - log2(y+1)
  return(l2r[isexpr])
}

library(ggplot2)
tpm_diff <- sapply(names(res_cntr),log2r)
df_diff <- melt(tpm_diff)

headtail <- function(x,n) c(head(x,n),tail(x,n))

ETP_DN_names <- names(headtail(sort(tpm_diff$DN_vs_ETP),10))
DPearly_vs_DN_names <- names(headtail(sort(tpm_diff$DPearly_vs_DN),10))
DPlate_vs_DPearly_names <- names(headtail(sort(tpm_diff$DPlate_vs_DPearly),10))
CD4SP_vs_DPlate_names <- names(headtail(sort(tpm_diff$CD4SP_vs_DPlate),10))
CD8SP_vs_DPlate_names <- names(headtail(sort(tpm_diff$CD8SP_vs_DPlate),10))
CD4SP_vs_CD8SP_names <- names(headtail(sort(tpm_diff$CD4SP_vs_CD8SP),10))

# plot boxplots
df_ETP_vs_DN <- subset(dplyr::select(data.frame(tpm_avg), 3, 6), rownames(tpm_avg) %in% ETP_DN_names)
df_DPearly_vs_DN <- subset(dplyr::select(data.frame(tpm_avg), 3, 4), rownames(tpm_avg) %in% DPearly_vs_DN_names)
df_DPlate_vs_DPearly <- subset(dplyr::select(data.frame(tpm_avg), 4, 5), rownames(tpm_avg) %in% DPlate_vs_DPearly_names)
df_CD4SP_vs_DPlate <- subset(dplyr::select(data.frame(tpm_avg), 5, 1), rownames(tpm_avg) %in% CD4SP_vs_DPlate_names)
df_CD8SP_vs_DPlate <- subset(dplyr::select(data.frame(tpm_avg), 5, 2), rownames(tpm_avg) %in% CD8SP_vs_DPlate_names)
df_CD4SP_vs_CD8SP <- subset(dplyr::select(data.frame(tpm_avg), 1, 2), rownames(tpm_avg) %in% CD4SP_vs_CD8SP_names)

df_ETP_vs_DN$gene <- rownames(df_ETP_vs_DN)
df_DPearly_vs_DN$gene <- rownames(df_DPearly_vs_DN)
df_DPlate_vs_DPearly$gene <- rownames(df_DPlate_vs_DPearly)
df_CD4SP_vs_DPlate$gene <- rownames(df_CD4SP_vs_DPlate)
df_CD8SP_vs_DPlate$gene <- rownames(df_CD8SP_vs_DPlate)
df_CD4SP_vs_CD8SP$gene <- rownames(df_CD4SP_vs_CD8SP)

df_ETP_vs_DN$direction <- ifelse(df_ETP_vs_DN$gene %in% names(head(sort(tpm_diff$DN_vs_ETP), 10)), yes = "Up", no = "Down")
df_DPearly_vs_DN$direction <- ifelse(df_DPearly_vs_DN$gene %in% names(head(sort(tpm_diff$DPearly_vs_DN), 10)), yes = "Up", no = "Down")
df_DPlate_vs_DPearly$direction <- ifelse(df_DPlate_vs_DPearly$gene %in% names(head(sort(tpm_diff$DPlate_vs_DPearly), 10)), yes = "Up", no = "Down")
df_CD4SP_vs_DPlate$direction <- ifelse(df_CD4SP_vs_DPlate$gene %in%  names(head(sort(tpm_diff$CD4SP_vs_DPlate), 10)), yes = "Up", no = "Down")
df_CD8SP_vs_DPlate$direction <- ifelse(df_CD8SP_vs_DPlate$gene %in%  names(head(sort(tpm_diff$CD8SP_vs_DPlate), 10)), yes = "Up", no = "Down")
df_CD4SP_vs_CD8SP$direction <- ifelse(df_CD4SP_vs_CD8SP$gene %in%  names(tail(sort(tpm_diff$CD4SP_vs_CD8SP), 10)), yes = "Up", no = "Down")

df_ETP_vs_DN$comparison <- "ETP_vs_DN"
df_DPearly_vs_DN$comparison <- "DN_vs_DPearly"
df_DPlate_vs_DPearly$comparison <- "DPearly_vs_DPlate"
df_CD4SP_vs_DPlate$comparison <- "DPlate_vs_CD4SP"
df_CD8SP_vs_DPlate$comparison <- "DPlate_vs_CD8SP"
df_CD4SP_vs_CD8SP$comparison <- "CD4SP_vs_CD8SP"

smash <- rbind(melt(df_ETP_vs_DN), melt(df_DPearly_vs_DN),melt(df_DPlate_vs_DPearly),melt(df_CD4SP_vs_DPlate),melt(df_CD8SP_vs_DPlate),melt(df_CD4SP_vs_CD8SP))


library(sleuth)
message("loading sleuth object")
so_thymo_transcripts <- sleuth_load("/media/colne37/hippopotamus/thymodevel/data/Turner/combined_transcripts_Turner_Thymo_sans_chrY_PAR_NM_AND_NR_ONLY_log2_sleuth")

# Collapse to gene
message("lollapse to gene")
refseq <- unique(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/hg38_annotation_for_sleuth_sans_chrY_PAR.tsv",header=F))
txi_thymo <- summarizeSleuthToGene(so_thymo_transcripts,refseq)

Turner_tpm <- data.frame(txi_thymo$abundance)

df_Turner <- dplyr::select(Turner_tpm, "Turner_DN","Turner_DPE","Turner_DPL","Turner_ETP","Turner_SP4", "Turner_SP8")
df_Turner$gene <- rownames(df_Turner)
df_Turner <- melt(df_Turner)
df_Turner$cell <- gsub(df_Turner$variable, pattern = "Turner_", replacement = "")
df_Turner[df_Turner$cell == "DPE",]$cell <- "DPearly"
df_Turner[df_Turner$cell == "DPL",]$cell <- "DPlate"
df_Turner[df_Turner$cell == "SP4",]$cell <- "CD4SP"
df_Turner[df_Turner$cell == "SP8",]$cell <- "CD8SP"

df_Turner$disease <- "Turner"
df_Turner <- dplyr::select(df_Turner, -variable)


asdf <- merge(df_Turner, dplyr::select(smash, -value), by.x = c("gene", "cell"), by.y = c("gene", "variable"))


smash$disease <- "Normal"



colnames(smash) <- c("gene", "direction","comparison","cell", "TPM", "disease")
colnames(asdf) <- c("gene", "cell", "TPM", "disease", "direction", "comparison")

smash_asdf <- rbind(smash, asdf)

smash_asdf$direction_corrected <- "KEK"
smash_asdf[smash_asdf$direction == "Down",]$direction_corrected <- "Up"
smash_asdf[smash_asdf$direction == "Up",]$direction_corrected <- "Down"

library(rstatix)
library(ggpubr)

stat.test <- compare_means(TPM ~ disease, smash_asdf, group.by = c("cell", "direction_corrected", "comparison"), method = "t.test", p.adjust.method = "BH")

stat.test <- stat.test %>%
  mutate(y.position = 400)

comparions_to_include <- c("DN_vs_DPearly", "DPlate_vs_CD4SP", "DPlate_vs_CD8SP", "CD4SP_vs_CD8SP")

stat.test_kek <- subset(dplyr::select(stat.test, -4, -p, -p.format, -p.signif, -method), comparison %in% comparions_to_include)

plot_df <- smash_asdf[smash_asdf$comparison %in% comparions_to_include,]

library(ggrepel)

filtzz <- intersect(unique(c(names(headtail(sort(tpm_diff$DPearly_vs_DN),10)),
                             names(headtail(sort(tpm_diff$CD4SP_vs_DPlate),10)),
                             names(headtail(sort(tpm_diff$CD8SP_vs_DPlate),10)),
                             names(headtail(sort(tpm_diff$CD4SP_vs_CD8SP),10))
)
), rownames(txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo$Celltype,function(y) mean(y>=1) )==1) ),]
))

genez_to_show <- c("TARP", "JCHAIN", "CD8A", "RORC", "CD8B", "RAG1", "RAG2", "CD40LG", "KLF2", "NOTCH3", "CCR7", "KLRK1", "CD4", "ZBTB7B", "HLA-F")

library(cowplot)
library(ggsci)

ggsave2("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_5B_left.pdf",width = 14,height = 8,
       ggplot(plot_df, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP")),
                           y=TPM)) + 
         geom_point(aes(col = disease)) + 
         geom_line(aes(group=interaction(gene, disease), col = disease)) + 
         facet_wrap(factor(direction_corrected, levels = c("Up", "Down"))~factor(comparison, levels = c("ETP_vs_DN","DN_vs_DPearly","DPearly_vs_DPlate","DPlate_vs_CD4SP","DPlate_vs_CD8SP","CD4SP_vs_CD8SP")), scales = "free", ncol = 4, nrow =2) + 
         scale_color_aaas() + 
         theme_AL_box() + 
         theme(legend.title = element_blank()) + 
         labs(x="", title = "Top 10 checkpoint genes per contrast")+
         geom_text_repel(
           data    = plot_df[plot_df$gene %in% genez_to_show,],
           mapping = aes(label = gene), min.segment.length = 0.11111111111111, nudge_y = 1, max.overlaps = 100)
)



###################
## Plot heatmaps ##
###################

source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")
source("/media/colne37/hippopotamus/thymodevel/my_scripts/summarizeSleuthToGene.R")
source("/media/colne37/hippopotamus/thymodevel/my_scripts/get_thymocyte_time_distance_from_pca.R")

#Import metadata
library(data.table)
message("loading metadata")
meta_turner <- fread("/media/colne37/hippopotamus/thymodevel/data/Turner/Turner_meta.tsv")
meta_turner <- meta_turner[order(meta_turner$Sample),]
message("final fixes")
meta_turner$Pseudotime <- 100
meta_turner[meta_turner$Celltype == "ETP"]$Pseudotime <- 1.0
meta_turner[meta_turner$Celltype == "DN"]$Pseudotime <- 2.0
meta_turner[meta_turner$Celltype == "DPE"]$Pseudotime <- 3.0
meta_turner[meta_turner$Celltype == "DPL"]$Pseudotime <- 4.0
meta_turner[meta_turner$Celltype == "SP4"]$Pseudotime <- 4.9
meta_turner[meta_turner$Celltype == "SP8"]$Pseudotime <- 5.1

meta_turner[meta_turner$Celltype == "DPE"]$Celltype <- "DPearly"
meta_turner[meta_turner$Celltype == "DPL"]$Celltype <- "DPlate"
meta_turner[meta_turner$Celltype == "SP4"]$Celltype <- "CD4SP"
meta_turner[meta_turner$Celltype == "SP8"]$Celltype <- "CD8SP"
meta_turner <- dplyr::select(meta_turner, Sample,Celltype,Replicate,Pseudotime,Batch,Sex)

message("loading metadata")
meta_thymo <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/meta_RNAseq.tsv")
meta_thymo <- meta_thymo[order(meta_thymo$Sample),]
meta_thymo <- dplyr::select(meta_thymo, -ID, -Diagnosis, -Age)

meta_thymo <- rbind(meta_turner, meta_thymo)

#Load sleuth object for differential expression
library(sleuth)
message("loading sleuth object")
so_thymo_transcripts <- sleuth_load("/media/colne37/hippopotamus/thymodevel/data/Turner/combined_transcripts_Turner_Thymo_sans_chrY_PAR_NM_AND_NR_ONLY_log2_sleuth")

# Collapse to gene
message("lollapse to gene")
refseq <- unique(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/hg38_annotation_for_sleuth_sans_chrY_PAR.tsv",header=F))
txi_thymo <- summarizeSleuthToGene(so_thymo_transcripts,refseq)

# Add metadata to sleuth object
so_thymo_transcripts$sample_to_covariates <- cbind(so_thymo_transcripts$sample_to_covariates,meta_thymo[,-1])


plot.heat.cntr <- function(df,res_nm,row_nm,cutoff=1e-3, ... ){
  require(ComplexHeatmap)
  require(RColorBrewer)
  require(circlize)
  strs <- strsplit(res_nm,"_")
  nm <- c(sapply(strs,"[[",1),sapply(strs,"[[",3))
  nm_reg <- paste0(nm,collapse="|")
  res <- res_cntr[[res_nm]]
  gn <- res$gene_id[res$qval < cutoff]
  idx <- intersect(gn,row.names(df))
  idx2 <- tpm_avg[idx,nm[1]]>=1|tpm_avg[idx,nm[2]]>=1
  df_red <- df[idx[idx2],grepl(nm_reg,colnames(df))]
  df_z <- t(apply(df_red,1,function(x) (x-mean(x))/sd(x) ))
  if(!missing(row_nm)) row.names(df_z)[!row.names(df_z) %in% row_nm] <- " "
  Heatmap(df_z, ... )
}

tpm_filter <- txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo$Celltype,function(y) mean(y>=1) )==1) ),]

colnames(tpm_filter)[colnames(tpm_filter) == 'Turner_DN'] <- 'TD2018_Turner_DN'
colnames(tpm_filter)[colnames(tpm_filter) == 'Turner_DPE'] <- 'TD2018_DPearly_Turner'
colnames(tpm_filter)[colnames(tpm_filter) == 'Turner_DPL'] <- 'TD2018_DPlate_Turner'
colnames(tpm_filter)[colnames(tpm_filter) == 'Turner_SP4'] <- 'TD2018_CD4SP_Turner'
colnames(tpm_filter)[colnames(tpm_filter) == 'Turner_SP8'] <- 'TD2018_CD8SP_Turner'


#tpm_avg <- t(apply(tpm_filter,1,function(x) tapply(x,meta_thymo$Celltype, mean) ))


colnames(tpm_filter)

library(ComplexHeatmap)

order_colz_DN_DPearly <- c("TD2018_DN_Rep1","TD2018_DN_Rep2","TD2018_DN_Rep3", "TD2018_DN_Rep4", "TD2018_DN_Rep5", "TD2018_Turner_DN","TD2018_DPearly_Rep0","TD2018_DPearly_Rep1","TD2018_DPearly_Rep2","TD2018_DPearly_Rep3", "TD2018_DPearly_Rep4", "TD2018_DPearly_Rep5", "TD2018_DPearly_Turner")

order_colz_DPearly_DPlate <- c("TD2018_DPearly_Rep0","TD2018_DPearly_Rep1","TD2018_DPearly_Rep2","TD2018_DPearly_Rep3", "TD2018_DPearly_Rep4","TD2018_DPearly_Rep5", "TD2018_DPearly_Turner",
                               "TD2018_DPlate_Rep0","TD2018_DPlate_Rep1","TD2018_DPlate_Rep2","TD2018_DPlate_Rep3", "TD2018_DPlate_Rep4","TD2018_DPlate_Rep5", "TD2018_DPlate_Turner")

order_colz_DPlate_CD4SP <- c("TD2018_DPlate_Rep0","TD2018_DPlate_Rep1","TD2018_DPlate_Rep2","TD2018_DPlate_Rep3", "TD2018_DPlate_Rep4","TD2018_DPlate_Rep5", "TD2018_DPlate_Turner",
                       "TD2018_CD4SP_Rep0","TD2018_CD4SP_Rep1","TD2018_CD4SP_Rep2","TD2018_CD4SP_Rep3", "TD2018_CD4SP_Rep4","TD2018_CD4SP_Rep5", "TD2018_CD4SP_Turner")

order_colz_DPlate_CD8SP <- c("TD2018_DPlate_Rep0","TD2018_DPlate_Rep1","TD2018_DPlate_Rep2","TD2018_DPlate_Rep3", "TD2018_DPlate_Rep4","TD2018_DPlate_Rep5", "TD2018_DPlate_Turner",
                             "TD2018_CD8SP_Rep0","TD2018_CD8SP_Rep1","TD2018_CD8SP_Rep2","TD2018_CD8SP_Rep3", "TD2018_CD8SP_Rep4","TD2018_CD8SP_Rep5", "TD2018_CD8SP_Turner")

order_colz_CD4SP_CD8SP <- c("TD2018_CD4SP_Rep0","TD2018_CD4SP_Rep1","TD2018_CD4SP_Rep2","TD2018_CD4SP_Rep3", "TD2018_CD4SP_Rep4","TD2018_CD4SP_Rep5", "TD2018_CD4SP_Turner",
                             "TD2018_CD8SP_Rep0","TD2018_CD8SP_Rep1","TD2018_CD8SP_Rep2","TD2018_CD8SP_Rep3", "TD2018_CD8SP_Rep4","TD2018_CD8SP_Rep5", "TD2018_CD8SP_Turner")

pdf("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_5B_right.pdf",width = 6,height = 8)

plot.heat.cntr(tpm_filter,"DPearly_vs_DN", genez_to_show, show_row_dend = F,col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="DPearly_vs_DN",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DN_DPearly)

plot.heat.cntr(tpm_filter,"CD4SP_vs_DPlate", genez_to_show,show_row_dend = F,col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD4SP_vs_DPlate",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DPlate_CD4SP)

plot.heat.cntr(tpm_filter,"CD8SP_vs_DPlate", genez_to_show,show_row_dend = F, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD8SP_vs_DPlate",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DPlate_CD8SP)

plot.heat.cntr(tpm_filter,"CD4SP_vs_CD8SP", genez_to_show,show_row_dend = F, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD4SP_vs_CD8SP",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_CD4SP_CD8SP)
dev.off()





############ Old way, below adding gene names to HM that are the strongest up & down in each contrast.

#pdf("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_5D_HM.pdf",width = 8,height = 12)

plot.heat.cntr(tpm_filter,"DPearly_vs_DN", names(headtail(sort(tpm_diff$DPearly_vs_DN),10)),show_row_dend = F,col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="DPearly_vs_DN",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DN_DPearly)

plot.heat.cntr(tpm_filter,"DPlate_vs_DPearly", names(headtail(sort(tpm_diff$DPlate_vs_DPearly),10)),show_row_dend = F, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="DPlate_vs_DPearly",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DPearly_DPlate)

plot.heat.cntr(tpm_filter,"CD4SP_vs_DPlate", names(headtail(sort(tpm_diff$CD4SP_vs_DPlate),10)),show_row_dend = F,col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD4SP_vs_DPlate",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DPlate_CD4SP)

plot.heat.cntr(tpm_filter,"CD8SP_vs_DPlate", names(headtail(sort(tpm_diff$CD8SP_vs_DPlate),10)),show_row_dend = F, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD8SP_vs_DPlate",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DPlate_CD8SP)

plot.heat.cntr(tpm_filter,"CD4SP_vs_CD8SP", names(headtail(sort(tpm_diff$CD4SP_vs_CD8SP),10)),show_row_dend = F, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD4SP_vs_CD8SP",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_CD4SP_CD8SP)
#dev.off()







