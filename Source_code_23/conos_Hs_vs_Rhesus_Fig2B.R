
# load("/domus/h1/dmiuso/private/Conos_HS/Conos_Hs_vs_Rhesus_INPUT_fig2B.RData")

# =============================================================

library(conos)
library(ggplot2)
library(pagoda2)
library(cowplot)
library(parallel)
#library(Matrix)
library(magrittr)
library(plyr)
library(tidyr)
library(dplyr)

new.label.info <- con.1$propagateLabels(labels = H.clust.TMP.from, method = mtd.name)

Lable.probab <- new.label.info[[3]]

www <- as.data.frame(as.matrix(Lable.probab))
www$cells.id <- rownames(www)
www <- www[(www$cells.id %in% names(H.clust.TMP.to)),]
meta.df <- data.frame("X" = names(H.clust.TMP.to), "clust" = as.character(H.clust.TMP.to), stringsAsFactors = F)
xxx <- www %>% pivot_longer(-cells.id, names_to = "Sharma.clust", values_to = "Probab")
XXX.long <- data.frame(xxx, "clust" = as.factor(mapvalues(xxx$cells.id, from = meta.df$X, to = meta.df$clust)))

pp <- ggplot(XXX.long, mapping = aes(x = Sharma.clust, y = Probab, fill = Sharma.clust)) +
  geom_violin(scale = "width", trim = T, show.legend = F) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  #theme.UD +
  facet_wrap(~ clust, ncol =1,
             #strip.position = "top",
             strip.position = "left",
             labeller = ) +
  theme(strip.text.y = element_text(size = 6, colour = "black", angle = 90))


pp <- pp + theme(plot.title = element_text(size = 10)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                                                              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pp