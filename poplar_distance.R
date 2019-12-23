library(reshape2)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(showtext)
library(purrr)
library(cowplot)
library(tidyr)
library(colorspace)
library(agricolae)

font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Md.otf",
         bolditalic = "/prop_fonts/prop_fonts/30. Minion     [1990 - Robert Slimbach]/MinionPro-BoldIt.otf")
showtext_auto()


#### remove previous statistic files ####
file.remove("stats_OD_poplarFW.csv")
file.remove("stats_hue_poplarFW.csv")


#### establish tukey-test functions ####
print.HSD.OD <- function(x) {
  aov1 <- aov(mean.od ~ bin + genotype, data = x)
  groups <- HSD.test(aov1, c("bin", "genotype"), alpha = 0.05)
  groups[["groups"]][["cell.type"]] <-  unique(x$cell.type)
  groups[["groups"]][["mean.od"]] <- -0.005
  write.table(
    groups[["groups"]],
    file = "stats_OD_poplarFW.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
}
print.HSD.hue <- function(x) {
  aov1 <- aov(mean.hue ~ bin + genotype, data = x)
  groups <- HSD.test(aov1, c("bin", "genotype"), alpha = 0.05)
  groups[["groups"]][["cell.type"]] <-  unique(x$cell.type)
  groups[["groups"]][["mean.hue"]] <- 290
  write.table(
    groups[["groups"]],
    file = "stats_hue_poplarFW.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
}

#### import poplar measurements ####
poplar <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/poplar_foodweb.csv")
poplar <- poplar[, -16]
poplar$replicate <- factor(poplar$replicate)
poplar$technical <- factor(poplar$technical)
poplar$number <- row(poplar)
poplar$cell.type <- plyr::revalue(poplar$cell.type,
                                  c(
                                    "F" = "Fibre",
                                    "V" = "Vessel",
                                    "R" = "Ray",
                                    "CB" = "Cambium"
                                  ))
poplar$adj.cell.type <- plyr::revalue(
  poplar$adj.cell.type ,
  c(
    "F" = "Fibre",
    "V" = "Vessel",
    "R" = "Ray",
    "CB" = "Cambium",
    "PA" = "Parenchyma"
  )
)


#### calculate distance to the cambium reference line 
# 5.9 is the number of pixels per µm
# X and Y are already measured according to scale, the ref values are given in pixels (see Fiji macro) ####
poplar$Distance <-
  apply(poplar[, c("X", "Y", "ref.x1", "ref.x2", "ref.y1", "ref.y2")],
        1 ,
        function(x) {
          a <- c(x[1], x[2])
          b <- c((x[3] / 5.9), (x[5] / 5.9))
          c <- c((x[4] / 5.9), (x[6] / 5.9))
          v1 <- b - c
          v2 <- a - b
          m <- cbind(v1, v2)
          d <- abs(det(m)) / sqrt(sum(v1 * v1))
          d
        })

poplar$diff <- poplar$OD.stained - poplar$OD.unstained
poplar.bg <-
  subset(
    poplar,
    cell.type == "Cambium" & adj.cell.type == "Cambium",
    select = c("genotype", "replicate", "technical", "diff")
  )
colnames(poplar.bg)[4] <- "OD.bg"

poplar.bg <- poplar.bg %>%
  group_by(genotype, replicate, technical) %>%
  mutate(OD.bg = mean(OD.bg, na.rm = TRUE))

poplar <-
  merge(poplar,
        unique(poplar.bg),
        by = c("genotype", "replicate", "technical"))

poplar$diff.adj <- poplar$diff - poplar$OD.bg
poplar$genotype <-
  ordered(poplar$genotype, levels = c("WT", "c4h", "ccr"))


#### calculate the correct hue on the 360 point circular scale ####
poplar$hue <- ((poplar$H.stained + 128) / 255 * 360)

poplar <- poplar %>%
  group_by(genotype, replicate, technical) %>%
  mutate(rel.dist = Distance / max(Distance, na.rm = TRUE))

poplar <- poplar %>%
  group_by(genotype, replicate, technical, cell.type, adj.cell.type) %>%
  mutate(rel.od = diff.adj / max(diff.adj, na.rm = TRUE))

poplar <-
  subset(
    poplar,
    cell.type != "Cambium" &
      adj.cell.type != "Parenchyma" & adj.cell.type != "Cambium"
  )
poplar$cell.type <- as.character(poplar$cell.type)
poplar$adj.cell.type <- as.character(poplar$adj.cell.type)
# poplar <- subset(poplar, cell.type == adj.cell.type) # uncomment to select only self anjacent cell walls
poplar$cell.type <- as.factor(as.character(poplar$cell.type))

poplar.bin <- poplar %>%
  mutate(bin = cut(
    Distance,
    breaks = c(-Inf, 50, 100, Inf),
    labels = c("I", "II", "III")
  ))

poplar.bin.pre <- poplar.bin %>%
  group_by(genotype, bin, cell.type, replicate) %>%
  summarise(
    mean.od = mean(diff.adj, na.rm = TRUE),
    sd.od = sd(diff.adj, na.rm = TRUE),
    mean.hue = mean(hue, na.rm = TRUE),
    sd.hue = sd(hue, na.rm = TRUE)
  )

poplar.bin.avg <- poplar.bin.pre %>%
  group_by(genotype, bin, cell.type) %>%
  summarise(
    mean.od = mean(mean.od, na.rm = TRUE),
    sd.od = sd(sd.od, na.rm = TRUE),
    mean.hue = mean(mean.hue, na.rm = TRUE),
    sd.hue = sd(sd.hue, na.rm = TRUE)
  )

poplar.bin$bin <-
  ordered(poplar.bin$bin, levels = c("I", "II", "III"))


#### set graph colours according to averaged measurements per genotype/cell type ####
barcols <- poplar.bin.avg[, c(1, 2, 3, 4, 6)]
colnames(barcols)[4] <- 'S'
colnames(barcols)[5] <- 'H'
barcols[, 4] <- barcols[, 4] * 2.5
barcols[, 4] <- ifelse(barcols[, 4] > 1, 1, barcols[, 4])
barcols[, 4] <- ifelse(barcols[, 4] < 0.1, 0.1, barcols[, 4])
barcols$V <- 0.95
barcols <- hex(HSV(data.matrix(barcols[, c(5, 4, 6)])))
poplar.bin$cell <-
  as.factor(paste(poplar.bin$genotype, poplar.bin$cell.type, poplar.bin$bin))
poplar.bin.avg$cell <-
  as.factor(paste(
    poplar.bin.avg$genotype,
    poplar.bin.avg$cell.type,
    poplar.bin.avg$bin
  ))
poplar.bin.pre$cell <-
  as.factor(paste(
    poplar.bin.pre$genotype,
    poplar.bin.pre$cell.type,
    poplar.bin.pre$bin
  ))
names(barcols) <- poplar.bin.avg$cell


#### set statistical letters for absorbance ####
poplar.bin.pre %>%
  group_by(cell.type) %>%
  do(data.frame(print.HSD.OD(.)))
letters.OD.monol <-
  read.csv("stats_OD_poplarFW.csv",
           header = FALSE)
colnames(letters.OD.monol) <-
  c("genobin", "mean.od", "group", "cell.type")
letters.OD.monol <-
  letters.OD.monol %>% separate(genobin, c("bin", "genotype"), sep = ":")
letters.OD.monol$bin <-
  factor(letters.OD.monol$bin, levels = c("I", "II", "III"))
letters.OD.monol$cell <-
  as.factor(
    paste(
      letters.OD.monol$genotype,
      letters.OD.monol$cell.type,
      letters.OD.monol$bin
    )
  )


#### plot absorbance ####
poplar.OD_dist <-
  ggplot(data = poplar.bin.pre, aes(
    y = mean.od,
    x = genotype,
    group = bin,
    fill = cell
  )) +
  annotate("rect",
                fill = "grey95",
                ymin = -Inf,
                ymax = Inf,
                xmin = 0.5,
                xmax = 1.5) +
  annotate("rect",
                fill = "grey95",
                ymin = -Inf,
                ymax = Inf,
                xmin = 2.5,
                xmax = 3.5) +
  geom_bar(
    stat = "identity",
    colour = NA,
    data = poplar.bin.avg,
    position = position_dodge(.9),
    alpha = 1,
    size = 0.1,
    width = 0.75
  ) +
  geom_point(
    # data = poplar.bin.pre,
    # fill = NA,
    position = position_jitterdodge(.5),
    # width = 0.1,
    alpha = 0.5,
    shape = 19,
    size = 1,
    stroke = 0
  ) +
  scale_y_continuous(
    limits = c(-0.005, 0.35),
    breaks = c(0, 0.1, 0.2, 0.3),
    labels = c(" 0.0", " 0.1", " 0.2", " 0.3"),
    expand = expand_scale(mult = c(0.15, 0))
  ) +
  scale_x_discrete(labels = c("WT",
                              expression(paste(italic("C4H"), "-RNAi")),
                              expression(paste(italic("CCR"), "-RNAi")))) +
  theme_minimal() +
  theme(
    text = element_text(size = 6, family = "Helvetica"),
    legend.position = "none",
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(
      size = 6,
      colour = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic",
      size = 6
    ),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  ) +
  labs(y = "Absorbance") +
  geom_text(
    position = position_dodge(0.9),
    data = letters.OD.monol,
    aes(label = group, y = mean.od, x = genotype),
    family = "Helvetica",
    angle = 90,
    hjust = 1,
    size = 6 / (14 / 5)
  ) +
  facet_grid(~ cell.type) +
  scale_fill_manual(values = barcols, guide = FALSE)

pdf("poplar_OD.pdf", width = 3.4, height = 1.5)
poplar.OD_dist
dev.off()


#### set statistical letters for hue ####
poplar.bin.pre %>%
  group_by(cell.type) %>%
  do(data.frame(print.HSD.hue(.)))
letters.hue.monol <-
  read.csv("stats_hue_poplarFW.csv",
           header = FALSE)
colnames(letters.hue.monol) <-
  c("genobin", "mean.hue", "group", "cell.type")
letters.hue.monol <-
  letters.hue.monol %>% separate(genobin, c("bin", "genotype"), sep = ":")
letters.hue.monol$bin <-
  factor(letters.hue.monol$bin, levels = c("I", "II", "III"))
letters.hue.monol$cell <-
  as.factor(
    paste(
      letters.hue.monol$genotype,
      letters.hue.monol$cell.type,
      letters.hue.monol$bin
    )
  )


#### plot hue ####
poplar.hue_dist <-
  ggplot(data = poplar.bin.pre, aes(
    y = mean.hue,
    x = genotype,
    colour = bin,
    fill = cell
  )) +
  annotate("rect",
           fill = "grey95",
           ymin = -Inf,
           ymax = Inf,
           xmin = 0.5,
           xmax = 1.5) +
  annotate("rect",
           fill = "grey95",
           ymin = -Inf,
           ymax = Inf,
           xmin = 2.5,
           xmax = 3.5) +
  annotate("text",
           x = c(0.65, 1.65, 2.65),
           y = 380,
           colour = "darkblue",
           label = "I",
           family = "Helvetica",
           fontface = 4,
           vjust = 0,
           size = 6 / (14 / 5)) +
  annotate("text",
           x = c(0.95, 1.95, 2.95),
           y = 380,
           colour = "darkblue",
           label = "II",
           family = "Helvetica",
           fontface = 4,
           vjust = 0,
           size = 6 / (14 / 5)) +
  annotate("text",
           x = c(1.3, 2.3, 3.3),
           y = 380,
           colour = "darkblue",
           label = "III",
           family = "Helvetica",
           fontface = 4,
           vjust = 0,
           size = 6 / (14 / 5)) +
  geom_jitter(
    data = poplar.bin.pre,
    position = position_jitterdodge(0.9),
    alpha = 0.75,
    shape = 21,
    size = 1,
    stroke = 0.25
  ) +
  stat_summary(
    fun.y = mean,
    fun.ymin = mean,
    fun.ymax = mean,
    geom = "crossbar",
    width = 0.5,
    size = 0.5,
    fatten = 0,
    position = position_dodge(0.9)
  ) +
  scale_y_continuous(
    limits = c(275, 390),
    breaks = c(300, 325, 350, 375),
    labels = c("300", "325", "350", "15")
  ) +
  scale_x_discrete(labels = c("WT",
                              expression(italic("c4h")),
                              expression(italic("ccr")))) +
  theme_minimal() +
  theme(
    text = element_text(size = 6, family = "Helvetica"),
    legend.position = "none",
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(
      size = 6,
      colour = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic",
      size = 6
    ),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  ) +
  labs(y = "Hue") +
  geom_text(
    position = position_dodge(0.9),
    data = letters.hue.monol,
    aes(label = group, y = mean.hue, x = genotype),
    family = "Helvetica",
    angle = 90,
    hjust = 1,
    size = 6 / (14 / 5)
  ) +
  scale_fill_manual(values = barcols, guide = FALSE) +
  scale_colour_manual(values = c("black", "black", "black")) +
  facet_grid(~ cell.type)

pdf("poplar_hue.pdf", width = 3.4, height = 1.5)
poplar.hue_dist
dev.off()

#### plot figure 7b,c ####
pdf("poplar_basic_grid.pdf", height = 4, width = 10)
plot_grid(
  poplar.OD_dist,
  poplar.hue_dist,
  labels = c('', ''),
  ncol = 1,
  nrow = 2,
  hjust = 0,
  vjust = 1,
  label_x = 0.003,
  label_y = 0.99,
  label_fontfamily = "Helvetica",
  rel_heights = c(1, 1)
)
dev.off()


####distance distribution and correlation with absorbance ####
pdf("distance_dist.pdf", width = 7, height = 4)
ggplot(data = poplar, aes(x = Distance, fill = genotype)) +
  geom_density(adjust = 2,
               color = NA,
               alpha = 0.75) +
  geom_vline(xintercept = 50,
             linetype = 2,
             color = "grey") +
  geom_vline(xintercept = 100,
             linetype = 2,
             color = "grey") +
  facet_grid(genotype ~ replicate) +
  scale_fill_few() +
  theme_few() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )
dev.off()

pdf("distance_od.pdf", width = 7, height = 7)
ggplot(data = subset(poplar, cell.type == "Vessel"),
       aes(x = Distance, y = rel.od, fill = genotype)) +
  geom_point(alpha = 0.75, shape = 21) +
  geom_smooth(method = "lm") +
  scale_x_continuous(limits = c(0, 300)) +
  facet_wrap( ~ genotype, ncol = 1) +
  scale_fill_few() +
  theme_few() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )
dev.off()