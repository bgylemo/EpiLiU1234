options(stringsAsFactors = F)

# Load data
source("/media/colne37/hippopotamus/thymodevel/my_scripts/data_prep_figure1/3_Thymocyte_load_sleuth_FC_NM_AND_NR_sans_chrY_PAR_transcripts.R")

# Filter out genes not expressed in any subtype
tpm_filter <- txi_thymo_transition$abundance[apply(txi_thymo_transition$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]

# Differential expression over timeseries
#so_thymo_transcripts <- sleuth_fit(so_thymo_transcripts, formula = ~Celltype,fit_name = "full")
#so_thymo_transcripts <- sleuth_fit(so_thymo_transcripts, formula =  ~1,fit_name = "reduced")
#so_thymo_transcripts <- sleuth_lrt(so_thymo_transcripts,"reduced","full")
#res_lrt <- sleuth_results(so_thymo_transcripts,'reduced:full',test_type = "lrt")

#write.table(res_lrt, file = "/media/god/my_book_8GB/LEGACY_DATA_BACKUP/TALL/res_lrt.tsv")

res_lrt <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/res_lrt.tsv")[,-1]

# Calculate Z-scores
df_z <- t(apply(tpm_filter[row.names(tpm_filter) %in% res_lrt$target_id,],1,function(x) (x-mean(x))/sd(x) ))

## Decide cluster n
#library(factoextra)
#library(cluster)
#pdf("plots/Thymocyte_timeseries_clustering_nclust_all.pdf")
#gridExtra::grid.arrange(
#  fviz_nbclust(df_z, kmeans, method = "wss",k.max = 25, nstart=15, iter.max=100),
#  fviz_nbclust(df_z, kmeans, method = "silhouette",k.max = 15, nstart=25,iter.max=100),
#  fviz_nbclust(df_z, kmeans, method = "gap_stat",nboot=100, k.max = 15, iter.max=100, nstart=25)
#)
#dev.off()

### NOT RUN: Test cluster n
library(ComplexHeatmap)
#Heatmap(df_z,km=4, show_row_names=F) 
#Heatmap(df_z,km=5, show_row_names=F)
#Heatmap(df_z,km=6, show_row_names=F)
#Heatmap(df_z,km=7, show_row_names=F)
#Heatmap(df_z,km=8, show_row_names=F)
#Heatmap(df_z,km=9, show_row_names=F)
#Heatmap(df_z,km=10, show_row_names=F)
#Heatmap(df_z,km=11, show_row_names=F)
#sapply(6:12,function(x){
#  cl <- kmeans(df_z,x,nstart = 25,iter.max = 100)$cluster
#  capture.output(split(names(cl),cl), file = paste0("OUT/GO/Clusters/",x,"/Thymocyte_timeseries_clusters.txt"))
#  sapply(seq_along(split(names(cl),cl)),function(i) write.table(split(names(cl),cl)[[i]],paste0("OUT/GO/Clusters/",x,"/Thymocyte_timeseries_cluster",i,".txt"),quote = F,sep = "\t",row.names = F,col.names = F) )
#  print(plot.cl(cl))
#})
### 7 clusters seems to yield the most information

# Cluster data into 7 clusters
set.seed(13) # for reproducible cluster order
cl <- kmeans(df_z,7,nstart = 25,iter.max = 100)$cluster

# split(names(cl),cl)

# Export list
#capture.output(split(names(cl),cl), file = "OUT/new/GO/Thymocyte_timeseries_clusters.txt")
#sapply(seq_along(split(names(cl),cl)),function(i) write.table(split(names(cl),cl)[[i]],paste0("OUT/new/GO/Thymocyte_timeseries_cluster",i,".txt"),quote = F,sep = "\t",row.names = F,col.names = F) )

df_z_for_HM <- df_z
colnamezz <- colnames(df_z_for_HM)
colnamezz_t <- gsub(colnamezz, pattern = "TD2018_", replacement = "")
colnames(df_z_for_HM) <- colnamezz_t

order_colz <- c("ETP_Rep1","ETP_Rep2","ETP_Rep3","ETP_Rep4","ETP_Rep5","DN_Rep1","DN_Rep2","DN_Rep3","DN_Rep4","DN_Rep5","DPearly_Rep0","DPearly_Rep1","DPearly_Rep2","DPearly_Rep3","DPearly_Rep4","DPearly_Rep5","DPlate_Rep0","DPlate_Rep1","DPlate_Rep2","DPlate_Rep3","DPlate_Rep4","DPlate_Rep5","CD4SP_Rep0","CD4SP_Rep1","CD4SP_Rep2","CD4SP_Rep3","CD4SP_Rep4","CD4SP_Rep5","CD8SP_Rep0","CD8SP_Rep1","CD8SP_Rep2","CD8SP_Rep3","CD8SP_Rep4","CD8SP_Rep5")

anno <-c(
"DNMT3B",
"HES1",
"HOXA9",
"LYL1",
"MYB",
"MYCN",
"NOTCH1",
"SPI1",
"TLX2",
"BCL11A",
"DNMT1",
"LMO2",
"RUNX1",
"STAT4",
"BACH2",
"CCR4",
"CCR5",
"CCR7",
"CD28",
"DNMT3A",
"ETS1",
"FOXP3",
"IKZF1",
"IL2RA",
"IRF4",
"MAF",
"NFATC1",
"NFATC2",
"RUNX3",
"STAT1",
"STAT2",
"STAT3",
"STAT5A",
"STAT5B",
"STAT6",
"TBX21",
"TNF",
"BCL2L1",
"NFATC3",
"RAG1",
"RORC",
"TCF12",
"BCL11B",
"CD3E",
"CD3G",
"E4F1",
"GATA3",
"KDM5A",
"LEF1",
"TET1")

#As we want to fix the order of the clusters, we have to re-order the gene-to-cluster assignment as a factor:
split <- factor(cl, levels=c("1","2","3","4","5","6","7"))

test_lol <- data.frame(namez = rownames(df_z_for_HM), index = 1:nrow(df_z_for_HM))

test_lol_kek <- test_lol[test_lol$namez %in% anno,]
  
ha = rowAnnotation(foo = anno_mark(at = test_lol_kek$index, labels = test_lol_kek$namez))

## Plot timeseries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
pdf("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_1D_right.pdf",width = 12,height = 16)
set.seed(13)
Heatmap(df_z_for_HM, split=split, show_row_names = F,column_order = order_colz, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))), show_row_dend = F, right_annotation = ha, use_raster = T, cluster_row_slices = FALSE, border = T)
dev.off()




plot.cl <- function(cl){
  ## Some ggplot madness ahead, basically plots the timeseries and splits SP stages into CD4/CD8
  require(ggplot2)
  require(reshape2)
  df_cl <- melt(df_z)
  df_cl$cluster <- cl
  df_cl$celltype <- sapply(strsplit(as.character(df_cl$Var2),"_"),"[[",2)
  df_cl$tp <- meta_thymo_transition$Pseudotime[match(df_cl$Var2,meta_thymo_transition$Sample)]
  ggplot() +
    stat_summary(data=df_cl[-grep("CD4",df_cl$celltype),],fun.data="mean_sdl",fun.args = list(mult=1),geom="ribbon", alpha=0.1,aes(y=value,x=round(tp),group=cluster,fill="CD8")) +
    stat_summary(data=df_cl[-grep("CD4",df_cl$celltype),],fun.y="mean",geom="line",alpha=0.5,aes(y=value,x=round(tp),group=cluster,col="CD8")) +
    stat_summary(data=df_cl[-grep("CD4",df_cl$celltype),],fun.y="mean",geom="point",alpha=0.5,aes(y=value,x=round(tp),group=cluster,col="CD8")) +
    stat_summary(data=df_cl[-grep("CD8",df_cl$celltype),],fun.data="mean_sdl",fun.args = list(mult=1),geom="ribbon", alpha=0.1,aes(y=value,x=round(tp),group=cluster, fill="CD4")) +
    stat_summary(data=df_cl[-grep("CD8",df_cl$celltype),],fun.y="mean",geom="line",alpha=0.5,aes(y=value,x=round(tp),group=cluster,col="CD4")) +
    stat_summary(data=df_cl[-grep("CD8",df_cl$celltype),],fun.y="mean",geom="point",alpha=0.5,aes(y=value,x=round(tp),group=cluster,col="CD4")) +
    facet_wrap(~cluster, ncol=1) +
    labs(y="Z-score") +
    scale_x_continuous(labels=c("1"="ETP\nCD34+\nCD1a-","2"="DN\nCD34+\nCD1a+","3"="DPearly\nCD3-","4"="DPlate\nCD3+","5"="SP\nCD3+")) +
    theme_AL_box(legend.position = "top", axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=4))+
    coord_cartesian(ylim=c(-2,2))+
    geom_hline(yintercept = 0, lty=2)
}

ggsave("/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_1D_left_CLUSTER.pdf",width = 6,height = 14, 
       plot_grid(rel_widths = c(0.75,1),plot.cl(cl),NULL))