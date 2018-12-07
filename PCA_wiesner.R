library(reshape2)
library(tidyr)
library(plyr)
library(dplyr)
library(factoextra)
library(ggfortify)
library(showtext)
library(cowplot)
library(ggthemes)

# import Helvetica Neue
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

###############################
# import measurements
###############################
phlog.monol <-
  read.csv(
    "/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
    skip = 2
  )
###############################
# calculate pixel values from OD
###############################
# phlog.monol$OD.stained <- 255/(10^phlog.monol$OD.stained)
# phlog.monol$OD.unstained <- 255/(10^phlog.monol$OD.unstained)

phlog.monol$genotype <-
  ordered(
    phlog.monol$genotype,
    levels = c(
      "col-0",
      "4cl1",
      "4cl2",
      "4cl1x2",
      "ccoaomt1",
      "fah1",
      "omt1",
      "ccr1-3",
      "ccr1xfah1",
      "cad4",
      "cad5",
      "cad4x5"
    )
  )

###############################
# set cell types according to measurement order
###############################
phlog.monol[1:50 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <- 
  "IF"
phlog.monol[51:100 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "MX"
phlog.monol[101:150 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "XF"
phlog.monol[151:200 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "PX"
phlog.monol[201:250 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "LP"
phlog.monol[251:300 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "PH"

###############################
# calculate the correct hue on the 360 point circular scale
###############################
phlog.monol$hue <- ((phlog.monol$h.stained + 128) / 255 * 360)

phlog.monol$replicate <-
  as.factor(as.character(phlog.monol$replicate))

###############################
# calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem
###############################
phlog.monol$diff <-
  phlog.monol$OD.stained - phlog.monol$OD.unstained
phlog.monol.bg <-
  ddply(
    subset(phlog.monol, cell.type == "PH", select = c(1, 2, 3, 4, 10)),
    c("genotype", "replicate", "technical"),
    summarise,
    OD.bg = mean(diff, na.rm = TRUE)
  )
phlog.monol.bg$cell.type <- NULL
phlog.monol <-
  merge(
    phlog.monol,
    phlog.monol.bg,
    all = TRUE,
    by = c("genotype", "replicate", "technical")
  )
phlog.monol$diff.adj <- phlog.monol$diff - phlog.monol$OD.bg
phlog.monol <- subset(phlog.monol, cell.type != "PH")

###############################
# average per replicate (for boxplots)
###############################
phlog.monol.pre <-
  ddply(
    phlog.monol,
    c("genotype", "cell.type", "replicate"),
    summarise,
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  )

###############################
# average per genotype (for barplots)
###############################
phlog.monol.avg <-
  ddply(
    phlog.monol.pre,
    c("genotype", "cell.type"),
    summarise,
    mean.hue2 = mean(mean.hue1, na.rm = TRUE),
    SD.hue2 = sd(mean.hue1, na.rm = TRUE),
    mean.OD2 = mean(mean.OD1, na.rm = TRUE),
    SD.OD2 = sd(mean.OD1, na.rm = TRUE)
  )


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

###############################
# shorten variable titles when only using OD
###############################
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

pca.wiesner.loading <- data.frame(pca.wiesner.post$rotation)
pca.wiesner.loading$genotype <- rownames(pca.wiesner.loading)
pca.wiesner.loading.melt <- melt(pca.wiesner.loading[, c(1,2,13)], id = "genotype")

pca.wiesner.loading.melt$genotype <-
  ordered(
    pca.wiesner.loading.melt$genotype,
    levels = c(
      "Col-0",
      "4cl1",
      "4cl2",
      "4cl1x4cl2",
      "ccoaomt1",
      "fah1",
      "omt1",
      "ccr1",
      "ccr1xfah1",
      "cad4",
      "cad5",
      "cad4xcad5"
    )
  )

pca.wiesner.ind <- data.frame(get_pca_ind(pca.wiesner.post)$cos2)
pca.wiesner.ind$cell.type <- celltypes.wiesner
pca.wiesner.ind <- tibble::rowid_to_column(pca.wiesner.ind, "number")
pca.wiesner.ind.melt <- melt(pca.wiesner.ind[, c(1,2,3,14)], id = c("number", "cell.type"))

gg.pca <- data.frame(pca.wiesner.post$x, cell.type = pca.wiesner$cell.type)
gg.rota <- data.frame(pca.wiesner.post$rotation)

pca <- ggplot(gg.pca, aes(x = PC1, y = PC2, fill = cell.type)) + 
  geom_hline(yintercept = 0, linetype = 1) +
  geom_vline(xintercept = 0, linetype = 1) +
  stat_ellipse(geom = "polygon", alpha = 0.5) +
  stat_ellipse(aes(fill = NULL), colour = "black", linetype = 2) +
  geom_point(shape = 21,
             size = 4,
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
    axis.title = element_text(size = 30),
    axis.text.y = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(
      size = 30,
      colour = "black",
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5),
    legend.position = c(0.1, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 30, colour = "black"),
    plot.margin = unit(c(2,2,2,2), "mm")
  ) 

pdf("PCA_wiesner.pdf")
pca
dev.off()

top <- plot_grid(
  fviz_pca_var(pca.wiesner.post,
               repel = TRUE,
               title = "",
               font.family = "Helvetica",
               col.var = "black",
               alpha.var = 0.4,
               ggtheme = theme_few()
  ),
  fviz_eig(pca.wiesner.post, 
           main = "", 
           font.family = "Helvetica",
           ggtheme = theme_few(),
           barcolor = NA,
           barfill = "#1d91c0"
           
  ),
  labels = c("(a)", "(b)"),
  label_fontfamily = "Helvetica",
  rel_widths = c(2,1),
  nrow = 1
)

mid <- plot_grid(
  fviz_contrib(pca.wiesner.post, 
               choice = "var", 
               axes = 1, 
               top = 12, 
               title = "", 
               font.family = "Helvetica",
               ggtheme = theme_few(),
               color = NA,
               fill = "#1d91c0"
  ),
  fviz_contrib(pca.wiesner.post, 
               choice = "var", 
               axes = 2, 
               top = 12, 
               title = "", 
               font.family = "Helvetica",
               ggtheme = theme_few(),
               color = NA,
               fill = "#1d91c0"
  ),
  labels = c("(c)", "(d)"),
  label_fontfamily = "Helvetica"
)

bottom <- ggplot(data = pca.wiesner.loading.melt, aes(x = genotype, y = value)) + 
  geom_point(shape = 21, size = 2) +
  geom_line(aes(group = variable, linetype = variable)) +
  scale_x_discrete(
    labels = c(
      "col-0",
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(paste(italic("4cl1"), "x", italic("4cl2"))),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(paste(italic("ccr1"), "x", italic("fah1"))),
      expression(italic("cad4")),
      expression(italic("cad5")),
      expression(paste(italic("cad4"), "x", italic("cad5")))
    )
  ) +
  labs(x = "Variable [genotype]",
       y = "Loadings") +
  theme_few() +
  theme(text = element_text(family = "Helvetica"),
        legend.title = element_blank(),
        legend.position = c(0.05,0.85))

pdf("PCA_supplemental.pdf", height = 13, width = 10)
plot_grid(
  top,
  mid,
  bottom,
  nrow = 3,
  labels = c("","","(e)"),
  label_fontfamily = "Helvetica",
  rel_heights = c(2,1,1)
)

# ggplot(data = pca.wiesner.ind.melt, aes(x = number, y = value)) + 
#   geom_point(aes(fill = cell.type, group = variable), shape = 21, size = 5) +
#   geom_line(aes(group = variable, linetype = variable))
dev.off()
