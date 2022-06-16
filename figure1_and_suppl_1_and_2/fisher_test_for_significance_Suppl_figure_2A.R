library(data.table)
library(dplyr)

sleuth <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/sans_chrY_PAR_NM_AND_NR_only_210427_all_beta_and_significance_transcripts.txt")

sleuth$significant <- ifelse(sleuth$pval < 0.00001 & (sleuth$b > 0.2 | sleuth$b < -0.2), yes = "sig", no = "ns")

# standard contigs
contigz <- c(paste0("chr", seq(1:22)), "chrX")

sleuth$type <- ifelse(sleuth$seqnames == "chrX", yes = "chrX", no = "autosome")

sleuth_count <- dplyr::select(sleuth[sleuth$seqnames %in% contigz,], target_id, significant, type) %>% dplyr::group_by( significant, type) %>% count()

sleuth_mat <- dcast(sleuth_count, type ~ significant, value.var = "n")
rownames(sleuth_mat) <- sleuth_mat$type
sleuth_mat$type <- NULL

test <- fisher.test(sleuth_mat)
test$p.value
