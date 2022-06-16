options(stringsAsFactors = F)
library(data.table)
library(immunarch)
library(cowplot)
source("/media/colne37/hippopotamus/thymodevel/my_scripts/AL_ggplot_themes.R")

# Read in TRUST4 data
immdata_sans_merged <- repLoad(.path = "/media/colne37/hippopotamus/thymodevel/data/Thymo/TCR_analysis/immunarch_sans_merged/")

down_samp <- repSample(
  immdata_sans_merged$data[13:40],
  .method = c("downsample"),
  .n = NA
)

ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Figure_4D.pdf", width = 4.5,
        plot_grid(repExplore(immdata_sans_merged$data, "volume") %>% vis(.test = F, .by = c("Cell", "Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6))))

ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Figure_4E.pdf", width = 2,
        plot_grid(repDiversity(down_samp, "d50") %>% vis(.test = F, .by = c("Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6)) + labs(x="")))

ggsave2(filename = "/media/colne37/hippopotamus/thymodevel/plots/Figure_4F.pdf", width = 4.5,
        plot_grid(  plot_grid(repClonality(down_samp, .method = "rare") %>% vis(.test = T, .by = c("Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6)))))