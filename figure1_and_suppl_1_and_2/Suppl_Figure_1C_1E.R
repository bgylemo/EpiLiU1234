#Load sleuth object
#source("/media/colne37/hippopotamus/thymodevel/my_scripts/data_prep_figure1/3_Thymocyte_load_sleuth_FC_NM_AND_NR_sans_chrY_PAR_transcripts.R")

##PLOT PCA
library(ggfortify)
plot.prcomp.loadings <- function(x,pc,n){
  require(ggplot2)
  lds <- sort(x$rotation[,pc],decreasing = T)
  lds_filt <- c(head(lds,n),tail(lds,n))
  qplot(x=reorder(names(lds_filt),-lds_filt,mean),y=lds_filt,geom="col") + labs(y="Contribution", title=paste("Principle component",pc)) + theme(axis.title.x = element_blank())
}

pca_thymo <- prcomp(t(txi_thymo_transition$counts))

# Read in stastistical test data
res_sig <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/res_lrt.tsv")[,-1]

# Filter out genes not expressed in any subtype
#tpm_filter <- txi_thymo_transition$abundance[apply(txi_thymo_transition$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]

#write.table(sort(pca_thymo$rotation[abs(pca_thymo$rotation[,1]) > 0.01,1],decreasing = T),"/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/PCA/thymocyte_PCA_component1.rnk",quote = F,sep = "\t",col.names = F)
#write.table(sort(pca_thymo$rotation[abs(pca_thymo$rotation[,2]) > 0.01,2],decreasing = T),"/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/PCA/thymocyte_PCA_component2.rnk",quote = F,sep = "\t",col.names = F)

#write.table(tpm_filter[row.names(tpm_filter) %in% res_sig$target_id[res_sig$qval < 1e-3],],"/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/thymocyte_timeseries.gct",sep = "\t",quote = F,col.names = NA)

#############################################
## Differential expression over timeseries ##
#############################################

## export for GSEA
export.gsea <- function(x,out){
  gsea_tab <- na.omit(x[,c('target_id','qval')])
  gsea_tab$qval <- -log10(gsea_tab$qval)
  write.table(gsea_tab,out,quote = F,sep = "\t",row.names = F,col.names = F)
}
#export.gsea(res_sig,"/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/thymocyte_timeseries.rnk")

#######################################
## Celltype-specific diff expression ##
#######################################

#so_thymo_transcripts <- sleuth_fit(so_thymo_transcripts, formula =  ~1,fit_name = "reduced")

#mod <- with(meta_thymo_transition,model.matrix(formula(~Celltype+0)))

#colnames(mod) <- gsub("Celltype","",colnames(mod))

sapply(colnames(mod),function(x){
  so_thymo_transcripts <<- sleuth_fit(so_thymo_transcripts,formula= ~mod[,x], fit_name=x)
  so_thymo_transcripts <<- sleuth_lrt(so_thymo_transcripts,"reduced",x)
})

#res_thymo_lrt_sig <- lapply(colnames(mod), function(x) sleuth_results(so_thymo_transcripts, paste0('reduced:',x), 'lrt', show_all = FALSE) )
#names(res_thymo_lrt_sig) <- colnames(mod)

## export relative expression
export.relexpr <- function(x,celltype,lrt,out,cutoff=1e-3){
  idx <- grepl(celltype,colnames(x))
  fc <- apply(x,1,function(y){
    df <- tapply(y,idx,mean)
    log2(df['TRUE']+1) - log2(df['FALSE']+1)
  })
  #filt <- subset(lrt,qval < cutoff)$target_id
  #fc_ord <- sort(fc[names(fc) %in% filt],decreasing = T)
  fc_ord <- sort(fc,decreasing = T)
  write.table(fc_ord,out,quote = F,sep = "\t",col.names = F)
}

#lapply(colnames(mod),function(x) export.relexpr(tpm_filter,x,res_thymo_lrt_sig[[x]],paste0("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/rel_expression/thymocyte_celltype_",x,".rnk")) )

# Export -log10(p) for GSEA
#lapply(names(res_thymo_lrt_sig),function(x) export.gsea(res_thymo_lrt_sig[[x]],paste0("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/pvals/thymocyte_celltype_",x,".rnk")) )



# Read in PCA GSEAPreranked output
library(data.table)
PCA1_neg_GSEA <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/GSEA_outputs/PCA_component_1_BP.GseaPreranked.1646223151067/gsea_report_for_na_neg_1646223151067.tsv")
PCA1_pos_GSEA <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/GSEA_outputs/PCA_component_1_BP.GseaPreranked.1646223151067/gsea_report_for_na_pos_1646223151067.tsv")

PCA1_neg_GSEA$NAME <- gsub(PCA1_neg_GSEA$NAME, pattern = "GOBP_", replacement = "")
PCA1_neg_GSEA$NAME <- gsub(PCA1_neg_GSEA$NAME, pattern = "_", replacement = " ")
PCA1_neg_GSEA$NAME <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", PCA1_neg_GSEA$NAME, perl=TRUE)
PCA1_pos_GSEA$NAME <- gsub(PCA1_pos_GSEA$NAME, pattern = "GOBP_", replacement = "")
PCA1_pos_GSEA$NAME <- gsub(PCA1_pos_GSEA$NAME, pattern = "_", replacement = " ")
PCA1_pos_GSEA$NAME <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", PCA1_pos_GSEA$NAME, perl=TRUE)


PCA2_neg_GSEA <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/GSEA_outputs/PCA_component_2_BP.GseaPreranked.1646223164829/gsea_report_for_na_neg_1646223164829.tsv")
PCA2_pos_GSEA <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GSEA/GSEA_outputs/PCA_component_2_BP.GseaPreranked.1646223164829/gsea_report_for_na_pos_1646223164829.tsv")

PCA2_neg_GSEA$NAME <- gsub(PCA2_neg_GSEA$NAME, pattern = "GOBP_", replacement = "")
PCA2_neg_GSEA$NAME <- gsub(PCA2_neg_GSEA$NAME, pattern = "_", replacement = " ")
PCA2_neg_GSEA$NAME <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", PCA2_neg_GSEA$NAME, perl=TRUE)
PCA2_pos_GSEA$NAME <- gsub(PCA2_pos_GSEA$NAME, pattern = "GOBP_", replacement = "")
PCA2_pos_GSEA$NAME <- gsub(PCA2_pos_GSEA$NAME, pattern = "_", replacement = " ")
PCA2_pos_GSEA$NAME <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", PCA2_pos_GSEA$NAME, perl=TRUE)

PCA1_components <- rbind(head(PCA1_neg_GSEA, 5), head(PCA1_pos_GSEA, 5))
PCA2_components <- rbind(head(PCA2_neg_GSEA, 5), head(PCA2_pos_GSEA, 5))

PCA1_order <- c(head(PCA1_neg_GSEA, 5)$NAME, head(PCA1_pos_GSEA, 5)$NAME)
PCA2_order <- c(head(PCA2_neg_GSEA, 5)$NAME, head(PCA2_pos_GSEA, 5)$NAME)


library(ggplot2)
library(cowplot)
source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")
ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_1C.pdf",
plot_grid(ncol=1,labels = c("PC1", "PC2"),
ggplot(PCA1_components, aes(x=NES,y=factor(NAME, levels = PCA1_order))) + geom_bar(stat="identity", width = 0.4) + coord_cartesian(xlim=c(-4,4)) + scale_x_continuous(position = "top", breaks = c(-4,-3,-2,-1,0,1,2,3,4)) + theme_AL_simple() + theme(axis.text.y = element_text(size=8)) + labs(y="", x="Normalised enrichment score (NES)") + geom_vline(xintercept = 0, lty=2),

plot_grid(ncol=2, NULL, rel_widths = c(0.125,1),
          ggplot(PCA2_components, aes(x=NES,y=factor(NAME, levels = PCA2_order))) + geom_bar(stat="identity", width = 0.4) + coord_cartesian(xlim=c(-3,3)) + scale_x_continuous(position = "top", breaks = c(-3,-2,-1,0,1,2,3)) + theme_AL_simple() + theme(axis.text.y = element_text(size=8)) + labs(y="", x="Normalised enrichment score (NES)") + geom_vline(xintercept = 0, lty=2))
))

cluster1 <- head(dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GO/PANTHER/cluster1.txt"), 1, 6), 10)
cluster2 <- head(dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GO/PANTHER/cluster2.txt"), 1, 6), 10)
cluster3 <- head(dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GO/PANTHER/cluster3.txt"), 1, 6), 10)
cluster4 <- head(dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GO/PANTHER/cluster4.txt"), 1, 6), 10)
cluster5 <- head(dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GO/PANTHER/cluster5.txt"), 1, 6), 10)
cluster6 <- head(dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GO/PANTHER/cluster6.txt"), 1, 6), 10)
cluster7 <- head(dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/new/GO/PANTHER/cluster7.txt"), 1, 6), 10)

cluster1$GO_short <- substring(cluster1$`GO biological process complete`, 1, nchar(cluster1$`GO biological process complete`)-13)
cluster2$GO_short <- substring(cluster2$`GO biological process complete`, 1, nchar(cluster2$`GO biological process complete`)-13)
cluster3$GO_short <- substring(cluster3$`GO biological process complete`, 1, nchar(cluster3$`GO biological process complete`)-13)
cluster4$GO_short <- substring(cluster4$`GO biological process complete`, 1, nchar(cluster4$`GO biological process complete`)-13)
cluster5$GO_short <- substring(cluster5$`GO biological process complete`, 1, nchar(cluster5$`GO biological process complete`)-13)
cluster6$GO_short <- substring(cluster6$`GO biological process complete`, 1, nchar(cluster6$`GO biological process complete`)-13)
cluster7$GO_short <- substring(cluster7$`GO biological process complete`, 1, nchar(cluster7$`GO biological process complete`)-13)

colnames(cluster1) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster2) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster3) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster4) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster5) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster6) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster7) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")

cluster1$Fold_enrichment <- as.numeric(cluster1$Fold_enrichment)
cluster1$GO_BP_short <- factor(cluster1$GO_BP_short, levels = rev(cluster1$GO_BP_short))
cluster1$cluster <- "Cluster 1"

cluster2$Fold_enrichment <- as.numeric(cluster2$Fold_enrichment)
cluster2$GO_BP_short <- factor(cluster2$GO_BP_short, levels = rev(cluster2$GO_BP_short))
cluster2$cluster <- "Cluster 2"

cluster3$Fold_enrichment <- as.numeric(cluster3$Fold_enrichment)
cluster3$GO_BP_short <- factor(cluster3$GO_BP_short, levels = rev(cluster3$GO_BP_short))
cluster3$cluster <- "Cluster 3"

cluster4$Fold_enrichment <- as.numeric(cluster4$Fold_enrichment)
cluster4$GO_BP_short <- factor(cluster4$GO_BP_short, levels = rev(cluster4$GO_BP_short))
cluster4$cluster <- "Cluster 4"

cluster5$Fold_enrichment <- as.numeric(cluster5$Fold_enrichment)
cluster5$GO_BP_short <- factor(cluster5$GO_BP_short, levels = rev(cluster5$GO_BP_short))
cluster5$cluster <- "Cluster 5"

cluster6$Fold_enrichment <- as.numeric(cluster6$Fold_enrichment)
cluster6$GO_BP_short <- factor(cluster6$GO_BP_short, levels = rev(cluster6$GO_BP_short))
cluster6$cluster <- "Cluster 6"

cluster7$Fold_enrichment <- as.numeric(cluster7$Fold_enrichment)
cluster7$GO_BP_short <- factor(cluster7$GO_BP_short, levels = rev(cluster7$GO_BP_short))
cluster7$cluster <- "Cluster 7"

smash_cluster <- rbind(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster7)

ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_1E.pdf", height = 20,
ggplot(smash_cluster, aes(x=Fold_enrichment,y=GO_BP_short)) + geom_bar(stat="identity", width = 0.4) +scale_x_continuous(position = "top", breaks = c(0,2,4,6,8,10,12,14,16)) + theme_AL_simple() + theme(axis.text.y = element_text(size=8)) + labs(y="", x="fold enrichment") + geom_vline(xintercept = 1, lty=2, col = "red") + coord_cartesian(xlim=c(0,16)) + facet_wrap(~cluster, ncol = 1, scale = "free")
)
