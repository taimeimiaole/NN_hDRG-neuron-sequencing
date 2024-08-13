
# load("/domus/h1/dmiuso/private/Conos_HS/Conos_Hs_vs_Sharma_INPUT_fig2A.RData")


mtd.name <- "diffusion"

new.label.info <- con.1$propagateLabels(labels = Sharma.meta.90.N, method = mtd.name)

Lable.probab <- new.label.info[[3]]

Name.extra <- "Raw.counts vs Sharma"

title.UD.probab <- paste(Name.extra, "\n", mtd.name, sep ="")

www <- as.data.frame(as.matrix(Lable.probab))

www$cells.id <- rownames(www)

www <- www[(www$cells.id %in% names(H.clust.TMP)),]

meta.df <- data.frame("X" = names(H.clust.TMP), "clust" = as.character(H.clust.TMP), stringsAsFactors = F)

xxx <- www %>% pivot_longer(-cells.id, names_to = "Sharma.clust", values_to = "Probab")

XXX.long <- data.frame(xxx, "clust" = as.factor(mapvalues(xxx$cells.id, from = meta.df$X, to = meta.df$clust)))

pp <- ggplot(XXX.long, mapping = aes(x = Sharma.clust, y = Probab, fill = Sharma.clust)) +
  geom_violin(scale = "width", trim = T, show.legend = F) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  facet_wrap(~ clust, ncol =1,
             #strip.position = "top",
             strip.position = "left",
             labeller = ) +
  theme(strip.text.y = element_text(size = 6, colour = "black", angle = 90))

pp + ggtitle(title.UD.probab) + theme(plot.title = element_text(size = 10))