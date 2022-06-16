options(stringsAsFactors = F)

# Load data
source("/media/colne37/hippopotamus/thymodevel/my_scripts/figure4_and_suppl_5/combined_Turner_Thymocyte_load_sleuth_FC_NM_AND_NR_sans_chrY_PAR_transcripts.R")

# Filter out genes not expressed in any subtype
tpm_filter <- txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo$Celltype,function(y) mean(y>=1) )==1) ),]

##PLOT PCA
library(ggfortify)
pca_thymo <- prcomp(t(txi_thymo$counts))

sample_order <- colnames(txi_thymo$counts)
meta_thymo$Sample <- factor(meta_thymo$Sample, levels = sample_order)

meta_thymo <- meta_thymo[order(meta_thymo$Sample),]
meta_thymo$Sex2 <- ifelse(meta_thymo$Sex == "T", yes = "Turner", no = "Other")

meta_thymo$Celltype <- factor(meta_thymo$Celltype, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP"))

PCA_plot <- autoplot(pca_thymo, data=meta_thymo, col = "Celltype", shape = "Sex", size=3, stroke=NA) + 
  theme_AL_box() + 
  annotate("text", x=0.2, y=-0.15, label= "ETP/DN")  +  
  annotate("text", x=-0.2, y=-0.05, label= "DPearly/DPlate")  +  
  annotate("text", x=0.1, y=0.15, label= "CD4SP/CD8SP")

library(cowplot)

ggsave2("/media/colne37/hippopotamus/thymodevel/plots/Figure_4B.pdf",width = 10,height = 8,
       plot_grid(PCA_plot))