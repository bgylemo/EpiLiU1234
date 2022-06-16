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

# Cleanup
rm(list=c("get.pca.dist"))