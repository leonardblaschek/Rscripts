library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
library(colorspace)
library(agricolae)
library(plyr)
library(dplyr)
library(rowr)
library(cowplot)

# import Helvetica Neue
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

# remove previous statistic files
file.remove("stats_OD_poplar.csv")
file.remove("stats_hue_poplar.csv")

# establish tukey-test functions
print.HSD.OD <- function(x) {
  aov1 <- aov(mean.OD1 ~ genotype, data = x)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups[["groups"]][["cell.type"]] <-  unique(x$cell.type)
  groups[["groups"]][["mean.OD1"]] <- -0.02
  write.table(
    groups[["groups"]],
    file = "stats_OD_poplar.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
}
print.HSD.hue <- function(x) {
  aov1 <- aov(mean.hue1 ~ genotype, data = x)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups[["groups"]][["cell.type"]] <-  unique(x$cell.type)
  groups[["groups"]][["mean.hue1"]] <- 308
  write.table(
    groups[["groups"]],
    file = "stats_hue_poplar.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
}

# calculate the correct hue on the 360 point circular scale
poplar$hue <- ((poplar$H.stained + 128) / 255 * 360)

poplar$replicate <-
  as.factor(as.character(poplar$replicate))

# average per replicate (for boxplots)
poplar.pre <-
  ddply(
    poplar,
    c("genotype", "cell.type", "replicate"),
    summarise,
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  )

# average per genotype (for barplots)
poplar.avg <-
  ddply(
    poplar.pre,
    c("genotype", "cell.type"),
    summarise,
    mean.hue2 = mean(mean.hue1, na.rm = TRUE),
    SD.hue2 = sd(mean.hue1, na.rm = TRUE),
    mean.OD2 = mean(mean.OD1, na.rm = TRUE),
    SD.OD2 = sd(mean.OD1, na.rm = TRUE)
  )

# set graph colours according to averaged measurements per genotype/cell type
barcols <- poplar.avg[, c(1, 2, 3, 5)]
colnames(barcols)[3] <- 'H'
colnames(barcols)[4] <- 'S'
barcols[, 4] <- barcols[, 4] * 2.5
barcols[, 4] <- ifelse(barcols[, 4] > 1, 1, barcols[, 4])
barcols[, 4] <- ifelse(barcols[, 4] < 0.1, 0.1, barcols[, 4])
barcols$V <- 0.95
barcols <- hex(HSV(data.matrix(barcols[, c(3, 4, 5)])))
poplar$cell <-
  as.factor(paste(poplar$genotype, poplar$cell.type))
poplar.avg$cell <-
  as.factor(paste(poplar.avg$genotype, poplar.avg$cell.type))
poplar.pre$cell <-
  as.factor(paste(poplar.pre$genotype, poplar.pre$cell.type))
names(barcols) <- poplar.avg$cell

# set statistical letters for hue
poplar.pre %>%
  group_by(cell.type) %>%
  do(data.frame(print.HSD.hue(.)))
letters.hue.monol <-
  read.csv("file:///home/leonard/R/Output/wiesner/stats_hue_poplar.csv",
           header = FALSE)
colnames(letters.hue.monol) <-
  c("genotype", "mean.hue1", "group", "cell.type")
letters.hue.monol$cell <-
  as.factor(paste(letters.hue.monol$genotype, letters.hue.monol$cell.type))

# plot hue
b <-
  ggplot(data = poplar.pre, aes(x = genotype, y = mean.hue1)) +
  geom_vline(xintercept = 1,
             size = 15,
             color = "grey95") +
  geom_vline(xintercept = 3,
             size = 15,
             color = "grey95") +
  # geom_jitter(
  #   data = poplar,
  #   aes(x = genotype, y = hue),
  #   width = 0.2,
  #   size = 1,
  #   stroke = 0.1,
  #   alpha = 0.3,
  #   shape = 21
  # ) +
  # geom_boxplot(
  #   outlier.size = 0,
#   position = position_dodge(width = 0.8),
#   size = 0.2,
#   width = 0.5
# ) +
geom_jitter(width = 0.1, alpha = 0.75, shape = 21, size = 3, stroke = 0.25) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.5, size = 0.5, fatten = 0) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    # axis.line.x = element_line(size = 0.75, lineend = "square"),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.y = element_line(size = 0.25, lineend = "square", color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      size =12,
      colour = "black",
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    legend.position = "none",
    strip.text = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  labs(y = "Hue", x = " ") +
  scale_x_discrete(
    labels = c(
      "WT",
      expression(italic("c4h")),
      expression(italic("ccr"))
    )
  ) +
  scale_y_continuous(
    breaks = c(310, 330, 350, 370),
    labels =
      c('310', '330', '350', '10'),
    limits = c(303, 370)
  ) +
  facet_wrap( ~ cell.type, ncol = 6) +
  aes(fill = cell) +
  scale_fill_manual(values = barcols) +
  geom_text(
    data = letters.hue.monol,
    aes(label = group),
    family = "Helvetica",
    angle = 0,
    colour = "black",
    hjust = 0.5,
    size = 4,
    position = position_dodge(width = 0.8)
  )
pdf("hue_monol_pop.pdf", height = 2.25, width = 5)
b
dev.off()

# set statistical letters for OD
poplar.pre %>%
  group_by(cell.type) %>%
  do(data.frame(print.HSD.OD(.)))
letters.OD.monol <-
  read.csv("file:///home/leonard/R/Output/wiesner/stats_OD_poplar.csv",
           header = FALSE)
colnames(letters.OD.monol) <-
  c("genotype", "mean.OD2", "group", "cell.type")
letters.OD.monol$cell <-
  as.factor(paste(letters.OD.monol$genotype, letters.OD.monol$cell.type))

# plot OD
a <- ggplot(poplar.avg, aes(x = genotype, y = mean.OD2)) +
  geom_vline(xintercept = 1,
             size = 15,
             color = "grey95") +
  geom_vline(xintercept = 3,
             size = 15,
             color = "grey95") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    legend.position = "none",
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.y = element_line(size = 0.25, lineend = "square", color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    ),
    plot.margin = unit(c(0,0,-0.4,0), "cm")
  ) +
  labs(y = "Absorbance", x = " ") +
  scale_x_discrete(breaks = c()) +
  scale_y_continuous(limits = c(-0.07, 0.55), breaks = c(0, 0.2, 0.4), labels = c(' 0.0', ' 0.2', ' 0.4')) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.85),
           width = 0.75) +
  facet_wrap( ~ cell.type, ncol = 6) +
  aes(fill = cell) +
  scale_fill_manual(values = barcols) + 
  # geom_violin(
  #   data = poplar,
  #   aes(x = genotype, y = diff.adj)
  # ) +
  # geom_boxplot(
  #   data = poplar.pre,
  #   aes(x = genotype, y = mean.OD1),
  #   outlier.size = 0.5,
  #   width = 0.4,
  #   position = position_dodge(width = 0.3),
  #   size = 0.2,
#   fill = "white"
# ) +
geom_jitter(data = poplar.pre,
            aes(x = genotype, y = mean.OD1), width = 0.1, alpha = 0.5) +
  geom_text(
    data = letters.OD.monol,
    aes(label = group),
    family = "Helvetica",
    angle = 0,
    colour = "black",
    hjust = 0.5,
    size = 4,
    position = position_dodge(width = 0.85)
  ) 
pdf("OD_monol_pop.pdf", height = 2, width = 5)
a
dev.off()

# plot grid
pdf("genotypes_poplar.pdf", height = 5, width = 5)
plot_grid(a,
          b,
          labels = c('A', 'B'),
          ncol=1,
          nrow = 2,
          hjust = 0,
          vjust = 1,
          label_x = 0.003,
          label_y = 0.99,
          label_fontfamily = "Helvetica",
          rel_heights = c(1, 1.1))
dev.off()