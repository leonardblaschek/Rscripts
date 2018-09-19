library(reshape2)
library(tidyr)
library(factoextra)
library(ggfortify)
library(showtext)
library(cowplot)

# import Helvetica Neue
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

# pca.wiesner <- phlog.monol.pre[, c(1:4,6)] # OD AND HUE
pca.wiesner <- phlog.monol.pre[, c(1:3,6)] # ONLY OD
pca.wiesner[1 + rep(seq(0, (nrow(pca.wiesner) - 1), by = 5), each = 1), 3] <- 
  1
pca.wiesner[2 + rep(seq(0, (nrow(pca.wiesner) - 1), by = 5), each = 1), 3] <- 
  2
pca.wiesner[3 + rep(seq(0, (nrow(pca.wiesner) - 1), by = 5), each = 1), 3] <- 
  3
pca.wiesner[4 + rep(seq(0, (nrow(pca.wiesner) - 1), by = 5), each = 1), 3] <- 
  4
pca.wiesner[5 + rep(seq(0, (nrow(pca.wiesner) - 1), by = 5), each = 1), 3] <- 
  5
pca.wiesner <- melt(pca.wiesner, id = c("cell.type", "replicate", "genotype"))
pca.wiesner$value[pca.wiesner$value < 0] <- 0.001
pca.wiesner <- dcast(pca.wiesner, cell.type + replicate ~ variable + genotype)

# shorten variable titles when only using OD
colnames(pca.wiesner) <- c("cell.type",
                           "replicate",
                           "Col-0",
                           "4cl1",
                           "4cl2",
                           "4cl1x4cl2",
                           "ccoaomt1",
                           "fah1",
                           "omt",
                           "ccr1",
                           "ccr1xfah1",
                           "cad4",
                           "cad5",
                           "cad4xcad5")


log.wiesner <- log(pca.wiesner[, -(1:2)])
celltypes.wiesner <- as.character(pca.wiesner[, 1])
pca.wiesner.post <- prcomp(log.wiesner, center = TRUE, scale. = TRUE)

gg.pca <- data.frame(pca.wiesner.post$x, cell.type = pca.wiesner$cell.type)
gg.rota <- data.frame(pca.wiesner.post$rotation)

pdf("PCA_wiesner.pdf")
p <- ggplot(gg.pca, aes(x = PC1, y = PC2, fill = cell.type)) + 
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.5) +
  geom_point(shape = 21,
             size = 3,
             stroke = 0.5,
             alpha = 0.75) +
  labs(x = "PC 1 (60%)",
       y = "PC 2 (26.1%)") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.ticks = element_line(size = 0.5, lineend = "square", color = "black"),
    axis.title = element_text(size = 20),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(
      size =20,
      colour = "black",
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5),
    legend.position = c(0.075, 0.875),
    legend.title = element_blank(),
    legend.text = element_text(size = 20, colour = "black"),
    plot.margin = unit(c(2,2,2,2), "mm")
  ) 
p
dev.off()

pdf("PCA_supplemental.pdf", height = 10, width = 10)
top <- plot_grid(
  fviz_pca_var(pca.wiesner.post,
               repel = TRUE,
               title = "",
               font.family = "Helvetica",
               col.var = "black",
               alpha.var = 0.4
  ),
  fviz_eig(pca.wiesner.post, main = "", font.family = "Helvetica"),
  labels = c("(a)", "(b)"),
  label_fontfamily = "Helvetica",
  rel_widths = c(2,1),
  nrow = 1
)

bottom <- plot_grid(
  fviz_contrib(pca.wiesner.post, choice = "var", axes = 1, top = 12, title = "", font.family = "Helvetica"),
  fviz_contrib(pca.wiesner.post, choice = "var", axes = 2, top = 12, title = "", font.family = "Helvetica"),
  labels = c("(c)", "(d)"),
  label_fontfamily = "Helvetica"
)
plot_grid(
  top,
  bottom,
  nrow = 2,
  labels = "",
  rel_heights = c(2,1)
)
dev.off()
print(pca.wiesner.post$rotation)
