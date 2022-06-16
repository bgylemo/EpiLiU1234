library(ggalluvial)
library(data.table)
library(dplyr)
source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")

# starting point
# read in the sleuth df
df_sleuth <- subset(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/sans_chrY_PAR_NM_AND_NR_only_210427_all_beta_and_significance_transcripts.txt"), seqnames %in% c("chrX"))

df_sleuth_short <- dplyr::select(df_sleuth, target_id, b, qval)

# load in our ASE data df, this is the starting point of the new classification. In this df I will add/change genes as escapees that Colm has eyeballed as low in methylation and that are sex-biased
df_ase <- subset(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/NM_NR_sans_chrY_PAR_ASE_table.tsv"), contig == "chrX" & sample == "F3" & keep == T)

# create classification df.
our_classification <- unique(dplyr::select(df_ase, gene, effectSize))
our_classification$category <- ifelse(our_classification$effectSize > 0.4, yes = "Inactive", no = "Escape")
our_classification <- unique(dplyr::select(our_classification, -effectSize))
our_classification <-our_classification[!grepl(our_classification$gene, pattern = ";"),]

ASE_classification <- subset(our_classification, !(gene == "PRKX" & category == "Inactive") & !(gene == "PUDP" & category == "Inactive") & !(gene == "TXLNG" & category == "Inactive"))

# add
smash <- merge(df_sleuth_short, ASE_classification, by.x="target_id", by.y = "gene", all.x = T)

# read in Colm's eyeball results
pois <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/Xlinked_promoter_methylation.csv", header = F)

pois$methylation_curated <- "Low"

smash1 <- merge(smash, pois, by.x = "target_id", by.y = "V1", all.x =T ) 

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

annotation_gene_table_short <- annotation_gene_table_short[!annotation_gene_table_short$name %in% c("NM_001320753", "NM_000351", "NM_001320754"),]

annotation_gene_table_short$start_fixed <- ifelse(annotation_gene_table_short$strand == "-", yes = annotation_gene_table_short$txEnd, no  = annotation_gene_table_short$txStart)

annotation_gene_table_short_tss1500 <- annotation_gene_table_short

anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)

anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)

anno_gr <- makeGRangesFromDataFrame(anno_short, ignore.strand=T, seqnames.field=c("CHR_hg38"), start.field=c("Start_hg38"),end.field=c("End_hg38"), keep.extra.columns = T)

annotation_gene_table_short_tss1500$first_tss1500 <- ifelse(annotation_gene_table_short_tss1500$strand == "+", yes = annotation_gene_table_short_tss1500$start_fixed - 500, no = annotation_gene_table_short_tss1500$start_fixed + 500)
annotation_gene_table_short_tss1500$second_tss1500 <- annotation_gene_table_short_tss1500$start_fixed + 1

annotation_gene_table_short_tss1500 <- transform(annotation_gene_table_short_tss1500, start_tss1500 = pmin(first_tss1500, second_tss1500))
annotation_gene_table_short_tss1500 <- transform(annotation_gene_table_short_tss1500, end_tss1500 = pmax(first_tss1500, second_tss1500))

annotation_gene_table_short_tss1500 <- dplyr::select(annotation_gene_table_short_tss1500, chrom, name2 ,name, start_tss1500, end_tss1500)

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

pois <- smash1$target_id


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

WOW_did_it_really_work[WOW_did_it_really_work$name2 == "SEPTIN6",]$name2 <- "SEPT6"

probes_methylation <- WOW_did_it_really_work[WOW_did_it_really_work$name2 %in% pois,]



melt_df_beta_pois <- merge(x=melt_df_beta, y=unique(probes_methylation[,-1]), by = c("probeID"))

melt_df_beta_pois_test <- melt_df_beta_pois %>% dplyr::group_by(name2) %>% dplyr::summarise(meanz = mean(value))

melt_df_beta_pois_test$meth_cat <- "00000"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz <= 0.25,]$meth_cat <- "0.25"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz > 0.25 & melt_df_beta_pois_test$meanz < 0.75,]$meth_cat <- "0.26-0.75"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz >= 0.75,]$meth_cat <- "0.76"


smash1$sexB <- ifelse(smash1$b >= 0.2 & smash1$qval < 0.05, yes = "sexB", no = "non_sexB")

smash2 <- merge(smash1, melt_df_beta_pois_test, by.x = "target_id", by.y = "name2", all.x=T)


df_ase <- subset(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/NM_NR_sans_chrY_PAR_ASE_table.tsv"), contig == "chrX" & sample == "F3" & keep == T)

# create classification df.
our_classification <- unique(dplyr::select(df_ase, gene, effectSize, cell))
our_classification_meanz <- our_classification %>% dplyr::group_by(gene) %>% dplyr::summarise(mean_ase = mean(effectSize))


smash2_short <- merge(dplyr::select(smash2[!is.na(smash2$meth_cat),], target_id, b, meth_cat),our_classification_meanz, by.x = "target_id", by.y = "gene" , all.x =T)

smash2_short[smash2_short$b > 1,]$b <- 1
smash2_short[smash2_short$b < -1,]$b <- -1

smash2_short_filt <- smash2_short

smash2_short_filt$stackz <- ifelse(smash2_short_filt$mean_ase > 0.4, yes = "Inactive", no ="Escape")

smash2_short_order <- smash2_short_filt[order(-meth_cat,b ), ]

PAR <- rbind(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/group-715_PAR1.csv"), fread("/media/colne37/hippopotamus/thymodevel/annotation_data/group-716_PAR2.csv"))
esc <- fread("/media/colne37/hippopotamus/thymodevel/annotation_data/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[esc$Gene_ID != "ENSG00000241489.3",]
esc <- esc[,-1]

chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"

chrx_annotation[chrx_annotation$Gene_name == "KAL1",]$Gene_name <- "ANOS1"
chrx_annotation[chrx_annotation$Gene_name == "CXorf36",]$Gene_name <- "DIPK2B"
chrx_annotation[chrx_annotation$Gene_name == "HDHD1",]$Gene_name <- "PUDP"
chrx_annotation[chrx_annotation$Gene_name == "RGAG4",]$Gene_name <- "RTL5"


smash2_short_melt <- melt(smash2_short_filt, id.vars = c("target_id", "stackz", "meth_cat"), measure.vars = c("b"))


ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Figure_3C.pdf", width = 18, height = 6,
        ggplot(smash2_short_melt[!(smash2_short_melt$target_id %in% PAR$`Approved symbol`),], aes(y=factor(meth_cat, levels = rev(c("0.25", "0.26-0.75","0.76"))), x=factor(target_id, levels = rev(smash2_short_order$target_id))) ) + 
          geom_tile(aes(fill=value)) + 
          scale_fill_gradientn(colours=c("red","grey","black"), limits =c(-1,1)) + 
          facet_wrap(stackz~., ncol = 1, scales = "free_x") + 
          theme_AL_box_rotX() + 
          theme(axis.text.x =element_text(size=4, face = "italic")) + labs(x="",y="")+
          geom_vline(xintercept = c("FIRRE", "UBA1", "SMARCA1", "PJA1", "TCEAL9", "FANCB", "NAP1L3", "ZRSR2", "RBBP7", "CACNA1F", "GCNA", "INE1", "MIR503HG")))


selectzzz <- smash2_short_melt[smash2_short_melt$stackz == "Inactive" | (is.na(smash2_short_melt$stackz) & smash2_short_melt$meth_cat != "0.25" & smash2_short_melt$value < 0.2),]

inactivezzz <- esc[esc$Gene_name %in% selectzzz$target_id,]
inactivezzz_fin <- inactivezzz[inactivezzz$Reported_XCI_status == "Escape",]


Figure3C_suppl_table <- smash2_short_melt
colnames(Figure3C_suppl_table) <- c("gene", "ASE_F3", "Methylation_category", "b", "log2FC") 

write.table(dplyr::select(Figure3C_suppl_table, -b), file = "/media/colne37/hippopotamus/thymodevel/plots/Suppl_table_9.tsv", quote = F, col.names = T, row.names = F, sep = "\t")



class_for_alluvium <- smash2_short_melt

class_for_alluvium[class_for_alluvium$meth_cat == "0.25" & class_for_alluvium$value >= 0.2 & is.na(stackz),]$stackz <- "Escape"
class_for_alluvium[class_for_alluvium$meth_cat != "0.25" & class_for_alluvium$value < 0.2 & is.na(stackz),]$stackz <- "Inactive"
class_for_alluvium[is.na(stackz),]$stackz <- "Unknown"








df_TPM <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/21_04_27_expression_matrix_unfiltered.tsv")
df_TPM_filt <- melt(df_TPM)

df_TPM_filt <- df_TPM_filt %>% dplyr::group_by(V1) %>% dplyr::summarise(meanz = mean(value))


# select the original Tuki df, this is the original classification of escapees
df_classification <- dplyr::select(chrx_annotation[chrx_annotation$Gene_name %in% class_for_alluvium$target_id,], Gene_name, category)
colnames(df_classification) <- c("gene", "category")
df_classification$axes <- "Tukiainen"


df_gylemo <- dplyr::select(class_for_alluvium, target_id, stackz)
colnames(df_gylemo) <- c("gene", "category")
df_gylemo$axes <- "Gylemo"


novelz <- data.frame(gene = setdiff(df_gylemo$gene, df_classification$gene),
                             category = "novel")
novelz$axes <- "Unannotated"


keklol <- df_gylemo[!df_gylemo$gene %in% df_classification$gene,]

df_novelz <- rbind(novelz, keklol)


is_lodes_form(df_novelz, key = "axes", value = "category", id = "gene")

df_novelz$axes <- factor(df_novelz$axes, levels = c("Unannotated", "Gylemo"))

df_novelz <- df_novelz[!df_novelz$gene %in% PAR$`Approved symbol` & df_novelz$gene %in% df_TPM_filt[df_TPM_filt$meanz > 1,]$V1,]


df_TPM_filt

smash_fin <- rbind(df_classification, df_gylemo)

smash_fin$axes <- factor(smash_fin$axes, levels = c("Tukiainen", "Gylemo"))

smash_fin <- smash_fin[!smash_fin$gene %in% PAR$`Approved symbol` & !smash_fin$gene %in% df_novelz$gene ,]

is_lodes_form(smash_fin, key = "axes", value = "category", id = "gene")


  


smashf_fin_filt <- smash_fin[smash_fin$gene %in% df_TPM_filt[df_TPM_filt$meanz > 1,]$V1,]

smashf_fin_filt$category <- factor(smashf_fin_filt$category, levels = c("Escape", "Inactive", "Variable", "Unknown"))

smashf_fin_filt %>% dplyr::group_by(category, axes) %>% dplyr::count()

smashf_fin_filt

figure3D_suppl_table <- smashf_fin_filt

figure3D_suppl_table[figure3D_suppl_table$axes == "Tukiainen",]$axes <- "previous assessment"
figure3D_suppl_table[figure3D_suppl_table$axes == "Gylemo",]$axes <- "thymocytes"

write.table(figure3D_suppl_table, file = "/media/colne37/hippopotamus/thymodevel/plots/Suppl_table_10.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Figure_3D.pdf",
        plot_grid(ggplot(smashf_fin_filt, aes(alluvium = gene, x = axes, stratum = category)) + 
                    geom_flow(aes(fill=category)) +
                    geom_stratum(aes(fill=category))  + 
                    geom_text(stat="count", aes(label=..count..))+
                    theme_AL_simple(legend.title = element_blank())+
                    scale_fill_manual(values = c("Escape" = "red", "Inactive" = "#869a9a", "Variable" = "purple","Unknown" = "blue"))
                  )
)



# how many expressed genes can we assess?
chrx_genes <- subset(fread("/media/colne37/hippopotamus/thymodevel/annotation_data/UCSC_refGene.tsv"), chrom %in% "chrX")

gois_covered <- subset(chrx_genes, name2 %in% df_TPM_filt[df_TPM_filt$meanz >= 1,]$V1)

unique(gois_covered$name2)

# excluding genes that does not have probes in their TSS (ending up with 362/521, covered/non-pseudo/expressed genes)
unique(gois_covered[gois_covered$name2 %in% df_gylemo[df_gylemo$category != "Unknown",]$gene,]$name2)



