################# Read in and prepare data ####################
# Load libraries
library(data.table)
library(pvclust)

# Read in data (created and saved above)
df_beta <- data.frame(fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/femme_only_methylation_beta_values_normalized.tsv", header = T))

# Read in metadata
metadata <- fread("/media/colne37/hippopotamus/thymodevel/data/res_and_expression_tables/thymocytes/meta_methylation.csv")

# melt
melt_df_beta <- melt(df_beta)

# calculate SD
variable_genes <- melt_df_beta %>% dplyr::group_by(probeID) %>% dplyr::summarise(sdz = sd(value))

# get the 1000 most variable genes
variable_genes_top1000 <- head(variable_genes[order(variable_genes$sdz, decreasing = T),], 1000)

df_pvclust <- df_beta[df_beta$probeID %in% variable_genes_top1000$probeID,]
rownames(df_pvclust) <- df_pvclust$probeID
df_pvclust$probeID <- NULL

# make pvclust object, including on the top 1000 most variable genes
results <- pvclust(df_pvclust, nboot = 100000, method = "euclidean", 
                   method.hclust = "complete")

pdf("/media/colne37/hippopotamus/thymodevel/plots/Figure_3A.pdf", width = 10)
plot(results)
pvrect(results, alpha=0.95)
dev.off()