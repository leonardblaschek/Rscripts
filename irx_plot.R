library(reshape2)
library(ggplot2)
library(ggthemes)
library(showtext)
library(dplyr)
library(cowplot)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

irx.plot <- function(x) {
irx.polygon <- ggplot(data = subset(irx.roi, genotype == x), aes(x = ROIx, y = ROIy, group = interaction(object, number))) +
  geom_polygon(colour = "black", aes(fill = Perim., linetype = object), alpha = 0.5, size = 0.5) +
  scale_fill_viridis_c(name = "Perimeter [µm]", limits = c(0, 100)) +
  guides(linetype = guide_legend(title="")) +
  xlim(0, 1.05) +
  ylim(0, 1.05) +
  labs(x = "Circularity", y = "Relative Distance") +
  theme(text = element_text(family = "Helvetica"),
        legend.position = c(0.025, 0.8),
        plot.margin = unit(c(0,0,0,0), "cm")) 

irx.polygon.densx <- ggplot(data = subset(irx.data, genotype == x), aes(Circ., linetype = object)) +
  geom_density() +
  xlim(0, 1.05) +
  # ylim(0, 1.05) +
  # theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank()) 

irx.polygon.densy <- ggplot(data = subset(irx.data, genotype == x), aes(Distance, linetype = object)) +
  geom_density() +
  xlim(0, 1.05) +
  # ylim(0, 1.05) +
  # theme_void() +
  theme(legend.position = "none",
        plot.margi = unit(c(0,0,1.25,0), "cm"),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank()) +
  coord_flip()

first.col <- plot_grid(irx.polygon, irx.polygon.densx, ncol = 1, rel_heights = c(10, 1), align = "v")
second.col <- plot_grid(irx.polygon.densy, NULL, ncol = 1, rel_heights = c(10, 1), align = "h")
scatter.grid <- plot_grid(first.col, second.col, ncol = 2, rel_widths = c(10, 1), align = "hv")
scatter.grid
}

pdf("irx_WT.pdf", width = 10, height = 10)
irx.plot("col-0")
dev.off()

pdf("irx_ccr1.pdf", width = 10, height = 10)
irx.plot("ccr1-3")
dev.off()

pdf("irx_ccoaomt1.pdf", width = 10, height = 10)
irx.plot("ccoaomt1")
dev.off()