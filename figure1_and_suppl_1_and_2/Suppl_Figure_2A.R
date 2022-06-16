source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")
options(stringsAsFactors = F)

library(data.table)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(karyoploteR)
library(rtracklayer)
library(cowplot)
library(ggplot2)
library(tidyr)


chromz <- c("chr1", "chr2", "chr3" ,"chr4", "chr5" , "chr6" ,"chr7" , "chr8" , "chr9" ,"chr10" , "chr11" , "chr12" , "chr13" ,"chr14" ,"chr15" ,"chr16" , "chr17", "chr18" , "chr19" , "chr20", "chr21" , "chr22" ,"chrX")

# Filter out genes not expressed
# Load data
TPM <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/21_04_27_expression_matrix_unfiltered.tsv")
tpm_melt <- melt(TPM)
tpm_melt$variable <- gsub(tpm_melt$variable, pattern = "TD2018_", replacement = "")
tpm_melt <- tpm_melt %>% separate(variable, into = c("cell", "sample"), sep = "_")

# Filter out genes not expressed in any subtype
tpm_filter <- tpm_melt %>% dplyr::group_by(V1, cell) %>% dplyr::summarise(meanz = mean(value))
tpm_filter <- subset(tpm_filter, meanz >= 1)


## Read in significance data
res_thymo_genes <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/sans_chrY_PAR_NM_AND_NR_only_210427_all_beta_and_significance_transcripts.txt")

# counts
nrow(subset(unique(dplyr::select(res_thymo_genes, target_id, pval)), pval < 0.00001))
nrow(subset(unique(dplyr::select(res_thymo_genes, target_id, pval, seqnames)), pval < 0.00001 & seqnames == "chrX"))

# calc -log10(p)
res_thymo_genes$log10pval <- -log10(res_thymo_genes$pval)

# Read in Tukiainen data frame
library(dplyr)
x.esc <- dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/landscape.Suppl.Table.13.csv"), Gene_name, Reported_XCI_status, Sex_bias_in_GTEx, XCI_across_tissues, XCI_in_single_cells)
esc <- dplyr::select(subset(x.esc), Gene_name, Reported_XCI_status)
esc <- esc[esc$Gene_name != ""]
esc[esc$Gene_name %in% "6-Sep"]$Gene_name <- "SEPT6"
esc <- esc[!(esc$Gene_name %in% "IDS" & esc$Reported_XCI_status == "Unknown")]


# merge our results with the Tukiainen df.
kary_res_thymo_genes <- merge(res_thymo_genes, esc, by.x = c("target_id"), by.y = c("Gene_name"), all.x = T)

# remove rows with NA contig, keep only columns we need.
kary_genes <- subset(dplyr::select(kary_res_thymo_genes, seqnames, start, end, strand, target_id, log10pval, Reported_XCI_status, qval, b), !is.na(seqnames) )

# Make all chrx genes not in Tukiainen paper unknown escape status.
kary_genes[is.na(kary_genes$Reported_XCI_status) & kary_genes$seqnames == "chrX",]$Reported_XCI_status <- "Unknown"

# I need to remove stuff that is not on the normal contigs
kary_genes <- subset(kary_genes, seqnames %in% chromz)

# Exclude NAs
kary_genes <- subset(kary_genes, !is.na(log10pval))

# extreme vals
extreme_vals <- subset(kary_genes, log10pval > 20)

# extreme vals
sign_genes <- subset(kary_genes, log10pval > 5)

# cap to 25
kary_genes[kary_genes$log10pval > 20]$log10pval <- 20

# make GRanges object
kary_genes_gr <- makeGRangesFromDataFrame(df = kary_genes, seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand", keep.extra.columns=TRUE)

# Remove chrY
kary_genes_gr_no_y <- kary_genes_gr[seqnames(kary_genes_gr) != "chrY"]

# check -log10(p) max values and round up to nearest 10 value.
library(plyr)
ymax_genes <- round_any(max(kary_genes_gr_no_y$log10pval), 5, f = ceiling)
#ymax_transcripts <- round_any(max(kary_transcripts_gr_no_y$log10pval), 10, f = ceiling)

# add ablines at -log10(5) and -log10(10). Calculate where to put it.
threshold_0.01_genes <- 10/ymax_genes
threshold_0.05_genes <- 5/ymax_genes

chromz <- c("chr1", "chr2", "chr3" ,"chr4", "chr5" , "chr6" ,"chr7" , "chr8" , "chr9" ,"chr10" , "chr11" , "chr12" , "chr13" ,"chr14" ,"chr15" ,"chr16" , "chr17", "chr18" , "chr19" , "chr20", "chr21" , "chr22" ,"chrX")

# Calculate how many genes are significant per chromosome (out of total genes we have for each chromosome).
library(dplyr)
kary_genes_count <- kary_genes[kary_genes$target_id %in% tpm_filter$V1,]
kary_genes_count$sig <- ifelse(kary_genes_count$log10pval > 5 & (kary_genes_count$b > 0.2 | kary_genes_count$b < -0.2), yes = "sig", no = "not_sig")
#kary_genes_count$sig <- ifelse(kary_genes_count$qval < 0.001, yes = "sig", no = "not_sig")
kary_genes_count_kek <- kary_genes_count %>% dplyr::group_by(seqnames, sig) %>% dplyr::count()
kary_genes_count_kek_wide <- dcast(kary_genes_count_kek, seqnames ~ sig, value.var="n")
#kary_genes_count_kek_wide[is.na(kary_genes_count_kek_wide$sig),]$sig <- 0
kary_genes_count_kek_wide$perc <- round(100*(kary_genes_count_kek_wide$sig / (kary_genes_count_kek_wide$not_sig + kary_genes_count_kek_wide$sig)), 2)
kary_genes_count_kek_wide$tot <- kary_genes_count_kek_wide$not_sig + kary_genes_count_kek_wide$sig

# read in chromosome lengths.
#chrom_length <- fread("/media/god/data1/R/thymodevel/chromosome_lengths.tsv")
#chrom_length$mid <- chrom_length$V2/2

#kary_genes_count_kek_wide <- merge(kary_genes_count_kek_wide, chrom_length, by.x = c("seqnames"), by.y = c("V1"))

DT.m1 = data.table::melt(kary_genes_count_kek_wide, id.vars = c("seqnames"))

library(ggsci)
library(scales)
#show_col(pal_npg("nrc")(10))


chrom_colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#7E6148FF","#B09C85FF","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#7E6148FF","#B09C85FF","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF", "#F39B7FFF")

chrom_colors_gg <- c("chr1" = "#E64B35FF","chr2" = "#4DBBD5FF","chr3" ="#00A087FF","chr4" ="#3C5488FF","chr5" ="#F39B7FFF","chr6" ="#8491B4FF","chr7" ="#91D1C2FF","chr8" ="#7E6148FF","chr9" ="#B09C85FF","chr10" ="#E64B35FF","chr11" ="#4DBBD5FF","chr12" ="#00A087FF","chr13" ="#3C5488FF","chr14" ="#F39B7FFF","chr15" ="#8491B4FF","chr16" ="#91D1C2FF","chr17" ="#7E6148FF","chr18" ="#B09C85FF","chr19" ="#E64B35FF","chr20" ="#4DBBD5FF","chr21" ="#00A087FF","chr22" ="#3C5488FF", "chrX" ="#F39B7FFF")

DT.m1$value <- round(DT.m1$value, 1)

chrom_colors_gg <- c("1" = "#E64B35FF","2" = "#4DBBD5FF","3" ="#00A087FF","4" ="#3C5488FF","5" ="#F39B7FFF","6" ="#8491B4FF","7" ="#91D1C2FF","8" ="#7E6148FF","9" ="#B09C85FF","10" ="#E64B35FF","11" ="#4DBBD5FF","12" ="#00A087FF","13" ="#3C5488FF","14" ="#F39B7FFF","15" ="#8491B4FF","16" ="#91D1C2FF","17" ="#7E6148FF","18" ="#B09C85FF","19" ="#E64B35FF","20" ="#4DBBD5FF","21" ="#00A087FF","22" ="#3C5488FF", "X" ="#F39B7FFF")

DT.m1$seqnames <- gsub(DT.m1$seqnames, pattern = "chr", replacement = "")

DT.m1$seqnames <- factor(DT.m1$seqnames, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"))


library(ggplotify)

pp <- getDefaultPlotParams(4)
pp$bottommargin <- 5
pp$rightmargin <- 0.0050
pp$leftmargin <- 0.04950

fig1b_manhattan <- as.ggplot(expression(kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL, chromosomes = chromz, plot.params = pp, labels.plotter = NULL),
                                        kpAxis(kp, ymin = 0, ymax = ymax_genes, tick.pos = c(0, 5, 10, 15, 20, 25)),
                                        kpPlotManhattan(kp, data = kary_genes_gr_no_y[kary_genes_gr_no_y$log10pval > 5], pval = kary_genes_gr_no_y[kary_genes_gr_no_y$log10pval > 5]$log10pval, logp = FALSE, ymin=0, ymax=ymax_genes, points.col = chrom_colors, suggestiveline = NULL, genomewideline = NULL, points.cex = 0.5),
                                        kpPlotManhattan(kp, data = kary_genes_gr_no_y[kary_genes_gr_no_y$log10pval < 5], pval = kary_genes_gr_no_y[kary_genes_gr_no_y$log10pval < 5]$log10pval, logp = FALSE, ymin=0, ymax=ymax_genes, points.col = chrom_colors, suggestiveline = NULL, genomewideline = NULL, points.cex = 0.5),
                                        kpAbline(kp, h = threshold_0.05_genes, col="black")))


fig1b_barplot <-
  ggplot(subset(DT.m1, variable == "perc" ), aes(x=seqnames,y=value,fill=seqnames, label = value)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=chrom_colors_gg) + 
  labs(y="", x = "") + 
  theme_AL_box() + 
  theme(legend.position = "none") + 
  geom_text(size = 2.5, vjust=0.15) + 
  coord_cartesian(ylim=c(0,5)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text.x = element_text(size=8), axis.text.y = element_text(size=8))



top_row <- plot_grid(fig1b_manhattan, fig1b_barplot, ncol = 1, rel_heights = c(1,0.5))


ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Suppl_Figure_2A.pdf", width = 150, 
        height = 100, units = "mm", 
        plot_grid(top_row)
)