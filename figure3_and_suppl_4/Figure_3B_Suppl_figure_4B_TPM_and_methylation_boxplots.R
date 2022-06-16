################# Read in and prepare data ####################
# Load libraries
library(data.table)
library(cowplot)
library(dplyr)
library(ggrepel)
library(rtracklayer)
library(GenomicRanges)
source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")

annotation_gene_table <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/UCSC_refGene.tsv")

chromz <- c(paste0("chr", 1:22), "chrX")

annotation_gene_table_short <- subset(dplyr::select(annotation_gene_table, chrom, txStart, txEnd, name2, name, strand), chrom %in% chromz)


annotation_gene_table_short$start_fixed <- ifelse(annotation_gene_table_short$strand == "-", yes = annotation_gene_table_short$txEnd, no  = annotation_gene_table_short$txStart)

annotation_gene_table_short_tss1500 <- annotation_gene_table_short



anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)

anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)

anno_gr <- makeGRangesFromDataFrame(anno_short, ignore.strand=T, seqnames.field=c("CHR_hg38"), start.field=c("Start_hg38"),end.field=c("End_hg38"), keep.extra.columns = T)



annotation_gene_table_short_tss1500$start_tss1500 <- annotation_gene_table_short_tss1500$start_fixed - 1500
annotation_gene_table_short_tss1500$end_tss1500 <- annotation_gene_table_short_tss1500$start_fixed + 500

annotation_gene_table_short_tss1500 <- annotation_gene_table_short_tss1500[,-c(2,3,6,7)]

colnames(annotation_gene_table_short_tss1500) <- c("seqnames","gene_id","transcript_id","start_tss1500","end_tss1500")

tss1500_gr <- makeGRangesFromDataFrame(annotation_gene_table_short_tss1500, seqnames.field=c("seqnames"), start.field=c("start_tss1500"), end.field=c("end_tss1500"), keep.extra.columns = T)

ol.f3 <- findOverlaps(anno_gr, tss1500_gr)
nm.f3 <- tapply(tss1500_gr$transcript_id[subjectHits(ol.f3)],queryHits(ol.f3),function(x) paste0(x,collapse=";") )

anno_short_tss1500 <- data.table(anno_short)
anno_short_tss1500[, transcript_id := NA]
anno_short_tss1500$transcript_id[as.numeric(names(nm.f3))] <- nm.f3

anno_short_tss1500 <- anno_short_tss1500[!is.na(anno_short_tss1500$transcript_id),]


testing_short <- dplyr::select(anno_short_tss1500, probeID, transcript_id)


luuul <- testing_short %>% tidyr::separate(transcript_id, sep=";" , into = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u"))


did_it_really_work <- dplyr::select(subset(melt(luuul, id.vars = "probeID"), !is.na(value)), -variable)


WOW_did_it_really_work <- unique(merge(dplyr::select(annotation_gene_table, name, name2), did_it_really_work, by.x = "name", by.y = "value", allow.cartesian=TRUE))

pois <- c("RAG1", "RORC", "CD3D", "CD3G", "CD8A")


df_beta <- data.frame(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/femme_only_methylation_beta_values_normalized.tsv", header = T))

# Read in metadata
metadata <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/meta_methylation.csv")

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



probes_methylation <- WOW_did_it_really_work[WOW_did_it_really_work$name2 %in% pois,]

# exclude RAG1 incorrect transcripts
probes_methylation <- probes_methylation[!probes_methylation$name %in% c("NM_001377277", "NM_001377278", "NM_001377279", "NM_001377280"),]

# exclude RORC incorrect transcript
probes_methylation <- probes_methylation[!probes_methylation$name %in% c("NM_005060"),]

# exclude CD8A incorrect transcripts
probes_methylation <- probes_methylation[!(probes_methylation$name2 == "CD8A" & !probes_methylation$name %in% c("NM_171827", "NM_001768", "NR_027353"))]


melt_df_beta_pois <- merge(x=melt_df_beta, y=unique(probes_methylation[,-1]), by = c("probeID"))


kek <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/21_04_27_expression_matrix_unfiltered.tsv")
kek_filt <- melt(kek[kek$V1 %in% pois,])
kek_filt <- kek_filt %>% tidyr::separate(variable, into = c(NA, "cell", "sample"), sep = "_")
kek_filt$sex <- ifelse(kek_filt$sample %in% c("Rep1", "Rep4", "Rep5"), yes = "female", no = "male")


tpmzzz <- 
  ggplot(kek_filt, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=value)) + 
  geom_boxplot() + 
  facet_wrap(~V1, scales = "free_y", ncol = 1) + 
  labs(x="", y="expression (TPM)") + 
  theme_AL_box_rotX() + 
  scale_y_continuous(limits = c(0, NA))


melt_df_beta_pois_median <- melt_df_beta_pois %>% dplyr::group_by(probeID, name2, cell) %>% dplyr::summarise(medianz = median(value))


methhh <- 
  ggplot(melt_df_beta_pois_median, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=medianz, col = probeID)) + 
  geom_point() +
  geom_line(aes(group=probeID)) +
  facet_wrap(~name2, ncol = 1) + 
  labs(x="", y="methylaion (median beta)") + 
  scale_y_continuous(limits = c(0, NA))+
  theme_AL_box_rotX(legend.position = "none")

ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Figure_3_meth_boxplots_and_TPMS.pdf", height = 14, width = 6,
plot_grid(ncol=2,
  tpmzzz, methhh)
)