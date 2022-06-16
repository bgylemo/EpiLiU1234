######### Figure 2C #############
library(dplyr)
library(data.table)

# read in Antonios plot themes #
source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")

cell_order <- c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP")

####### Read in data ########
ase.f3_unfilt <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/NM_NR_sans_chrY_PAR_ASE_table.tsv")

ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt[!is.na(ase.f3_unfilt$gene) & sample == "F3" & contig == "chrX"]

###### add PAR info #####
PAR <- rbind(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/group-715_PAR1.csv"), fread("/media/colne37/hippopotamus/thymodevel/annotation_data/group-716_PAR2.csv"))
PAR$GENES <- PAR$`Approved symbol`
ase.f3_unfilt_with_chrx_info <- merge(ase.f3_unfilt_with_chrx_info, PAR, by.x = "gene", by.y = "GENES", all.x = T)
ase.f3_unfilt_with_chrx_info[is.na(PAR)]$PAR <- "NON_PAR"
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$PAR == "PAR1" | ase.f3_unfilt_with_chrx_info$PAR == "PAR2"]$PAR <- "PAR"

##### add Esc info ######
x.esc <- dplyr::select(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/landscape.Suppl.Table.13.csv"), Gene_name, Reported_XCI_status, Sex_bias_in_GTEx, XCI_across_tissues, XCI_in_single_cells)
x.esc[x.esc$Gene_name == "6-Sep"]$Gene_name <- "SEPT6"
x.esc <- subset(x.esc, !(Reported_XCI_status == "Unknown" & Gene_name == "IDS"))
x.esc <- x.esc[x.esc$Gene_name != ""]

##### add info to df #########
ase.f3_unfilt_with_chrx_info <- merge(ase.f3_unfilt_with_chrx_info, x.esc, by.x = "gene", by.y = "Gene_name", all.x = T)

# remove overlapping hSNPs
ase.f3_unfilt <- ase.f3_unfilt[!grepl(ase.f3_unfilt$gene, pattern = ";")]
ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt_with_chrx_info[!grepl(ase.f3_unfilt_with_chrx_info$gene, pattern = ";")]

# change cell order
ase.f3_unfilt$cell <- factor(ase.f3_unfilt$cell, levels = cell_order)
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

################# Small fixes ####################
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Unknown"]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info[is.na(ase.f3_unfilt_with_chrx_info$Reported_XCI_status)]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info$status <- ifelse(ase.f3_unfilt_with_chrx_info$effectSize > 0.4, yes = "Monoallelic", no = "Biallelic")

ase.f3_unfilt$cell <- factor(ase.f3_unfilt$cell, levels = cell_order)
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

chrx_gene_order <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/chrx_gene_order.tsv")

## Idenfity 'escape' genes that does not escape in thymocytes
id_genez <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$keep == T & ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Escape",]

######### plot #########

library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)


point_size <- 1.5


####### Check unfiltered inactive genes throughout development ##########

unfilt_inac <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Inactive" & ase.f3_unfilt_with_chrx_info$gene != "SEPT6" & ase.f3_unfilt_with_chrx_info$variantCallDepth > 20 & ase.f3_unfilt_with_chrx_info$totalCount >= 10,]

length(unique(unfilt_inac$position))
length(unique(unfilt_inac$gene))


range(unfilt_inac$totalCount)
range(unfilt_inac$refCount)
range(unfilt_inac$altCount)

leeeel <- dplyr::select(unfilt_inac, altCount, refCount, gene, cell, position)

leeeel_min <- transform(leeeel, min = pmin(altCount, refCount))
leeeel_min_max <- transform(leeeel_min, max = pmax(altCount, refCount))

melt_leeeel_min_max <- melt(dplyr::select(leeeel_min_max, -altCount, -refCount), id.vars = c("gene", "cell", "position"), measure.vars = c("min", "max"))

melt_leeeel_min_max_countz <- melt_leeeel_min_max %>% dplyr::group_by(cell, variable) %>% dplyr::summarise(asdf = sum(value))

#inac_plot <- ggplot(melt_leeeel_min_max_countz, aes(x=cell, y=asdf, fill=variable, label = asdf)) + geom_bar(stat="identity", position = "dodge") + labs(title = "Inactive genes") + geom_text(position = position_dodge(width=.85))

stoopid <- as.data.frame(dcast(leeeel_min_max, gene+position~cell, value.var = c("min", "max")))

stoopid[is.na(stoopid)] <- 0

melt_leeeel_min_max_countz_dcast <- dcast(melt_leeeel_min_max_countz, cell ~ variable, value.var = "asdf")

melt_leeeel_min_max_countz_dcast$tot <- melt_leeeel_min_max_countz_dcast$min + melt_leeeel_min_max_countz_dcast$max

#######################







####### Check immune genes throughout development ##########

gois_1 <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/immune_genes_nri2815.tsv", header = F)
gois_2 <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/immune_genes_s13293-019-0278.tsv", header = T)
gois_3 <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/Immune_chrX_genes_fimmu.2017.01455.tsv", header = T)

genez <- unique(sort(c(gois_1$V1, gois_2$Gene_name, gois_3$gene)))

unfilt_immune <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$gene %in% genez & ase.f3_unfilt_with_chrx_info$variantCallDepth > 20 & ase.f3_unfilt_with_chrx_info$totalCount >= 10 & ase.f3_unfilt_with_chrx_info$effectSize > 0.4,]

ggplot(unfilt_immune, aes(x=cell, y = effectSize)) + facet_wrap(~gene) + geom_line(aes(group=gene)) + coord_cartesian(ylim=c(0,0.5))

range(unfilt_immune$totalCount)
range(unfilt_immune$refCount)
range(unfilt_immune$altCount)

leeeel2 <- dplyr::select(unfilt_immune, altCount, refCount, gene, cell, position)

leeeel_min2 <- transform(leeeel2, min = pmin(altCount, refCount))
leeeel_min_max2 <- transform(leeeel_min2, max = pmax(altCount, refCount))

melt_leeeel_min_max2 <- melt(dplyr::select(leeeel_min_max2, -altCount, -refCount), id.vars = c("gene", "cell", "position"), measure.vars = c("min", "max"))

melt_leeeel_min_max_countz2 <- melt_leeeel_min_max2 %>% dplyr::group_by(cell, variable) %>% dplyr::summarise(asdf = sum(value))


length(unique(melt_leeeel_min_max2$position))
length(unique(melt_leeeel_min_max2$gene))


stoopid2 <- as.data.frame(dcast(leeeel_min_max2, gene+position~cell, value.var = c("min", "max")))

stoopid2[is.na(stoopid2)] <- 0

melt_leeeel_min_max_countz_dcast2 <- dcast(melt_leeeel_min_max_countz2, cell ~ variable, value.var = "asdf")

melt_leeeel_min_max_countz_dcast2$tot <- melt_leeeel_min_max_countz_dcast2$min + melt_leeeel_min_max_countz_dcast2$max

# INCLUDED IMMUNE GENES
c(ARHGAP4,ARHGEF6,BTK,CD40LG,DOCK11,ELF4,G6PD,IL2RG,IL9R,ITM2A,MPP1,MSL3,PIGA,RPS6KA3,SASH3,TAB3,TFE3,TMSB4X,XIAP)

#######################
immune_plot <- ggplot(melt_leeeel_min_max_countz_dcast2, aes(x=cell, y=min, label = paste0("total: ", tot))) + geom_bar(stat="identity") + labs(x="",y="minor allele count", title = "87 immune genes; 32 hSNPs in 19 genes, wes>20, rnaseq>=10") + theme_AL_box() + geom_text(aes(y=10.5)) + geom_text(data=melt_leeeel_min_max_countz_dcast2, aes(label=min))

inac_plot <- ggplot(melt_leeeel_min_max_countz_dcast, aes(x=cell, y=min, label = paste0("total: ", tot))) + geom_bar(stat="identity") + labs(x="",y="minor allele count", title = "inactive (not SEPT6); 121 hSNPs in 68 genes, wes>20, rnaseq>=10") + theme_AL_box() + geom_text(aes(y=22.5)) + geom_text(data=melt_leeeel_min_max_countz_dcast, aes(label=min))

library(cowplot)
ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Figure_2G.pdf", height = 8, width=14,
        plot_grid(ncol=2, inac_plot, immune_plot))



