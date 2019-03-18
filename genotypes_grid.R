library(showtext)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
library(colorspace)
library(agricolae)
library(dplyr)
library(rowr)
library(cowplot)
library(png)
library(grid)


#### import Helvetica Neue ####
font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()


#### remove previous statistic files ####
file.remove("stats_OD.csv")
file.remove("stats_hue.csv")


#### establish tukey-test functions ####
print.HSD.OD <- function(x) {
  aov1 <- aov(mean.OD1 ~ genotype, data = x)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups[["groups"]][["cell.type"]] <-  unique(x$cell.type)
  groups[["groups"]][["mean.OD1"]] <- -0.007
  write.table(
    groups[["groups"]],
    file = "stats_OD.csv",
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
    file = "stats_hue.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
}


#### import measurements ####
phlog.monol <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
           skip = 2)

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


#### set cell types according to measurement order ####
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

phlog.monol$cell.type <- factor(phlog.monol$cell.type)

#### calculate the correct hue on the 360 point circular scale ####
phlog.monol$hue <- ((phlog.monol$h.stained + 128) / 255 * 360)

phlog.monol$replicate <-
  as.factor(as.character(phlog.monol$replicate))


#### calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem ####
phlog.monol$diff <-
  phlog.monol$OD.stained - phlog.monol$OD.unstained
phlog.monol.bg <- phlog.monol %>%
  filter(cell.type == "PH") %>%
  select(1:4, 10) %>%
  group_by(genotype, replicate, technical) %>%
  summarise(OD.bg = mean(diff, na.rm = TRUE))

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


#### average per replicate (for boxplots) ####
phlog.monol.pre <-  phlog.monol %>%
    group_by(genotype, cell.type, replicate) %>%
    summarise(
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  )


#### average per genotype (for barplots) ####
phlog.monol.avg <- phlog.monol.pre %>%
    group_by(genotype, cell.type) %>%
    summarise(
    mean.hue2 = mean(mean.hue1, na.rm = TRUE),
    SD.hue2 = sd(mean.hue1, na.rm = TRUE),
    mean.OD2 = mean(mean.OD1, na.rm = TRUE),
    SD.OD2 = sd(mean.OD1, na.rm = TRUE)
  )


#### set graph colours according to averaged measurements per genotype/cell type ####
barcols <- phlog.monol.avg[, c(1, 2, 3, 5)]
colnames(barcols)[3] <- 'H'
colnames(barcols)[4] <- 'S'
barcols <- barcols %>%
  mutate(S = S*2.5) %>%
  mutate(S = ifelse(S > 1, 1, S),
         S = ifelse(S < 0.1, 0.1, S))
barcols$V <- 0.95
barcols <- hex(HSV(data.matrix(barcols[, c(3, 4, 5)])))
phlog.monol$cell <-
  as.factor(paste(phlog.monol$genotype, phlog.monol$cell.type))
phlog.monol.avg$cell <-
  as.factor(paste(phlog.monol.avg$genotype, phlog.monol.avg$cell.type))
phlog.monol.pre$cell <-
  as.factor(paste(phlog.monol.pre$genotype, phlog.monol.pre$cell.type))
names(barcols) <- phlog.monol.avg$cell


#### set statistical letters for hue ####
phlog.monol.pre %>%
  group_by(cell.type) %>%
  do(data.frame(print.HSD.hue(.)))
letters.hue.monol <-
  read.csv("stats_hue.csv",
           header = FALSE)
colnames(letters.hue.monol) <-
  c("genotype", "mean.hue1", "group", "cell.type")
letters.hue.monol$cell <-
  as.factor(paste(letters.hue.monol$genotype, letters.hue.monol$cell.type))


#### plot hue ####
b <-
  ggplot(data = phlog.monol.pre, aes(x = genotype, y = mean.hue1)) +
  geom_vline(xintercept = 1,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 3,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 5,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 7,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 9,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 11,
             size = 5,
             color = "grey95") +
  geom_jitter(
    width = 0.1,
    alpha = 0.75,
    shape = 21,
    size = 3,
    stroke = 0.25
  ) +
  stat_summary(
    fun.y = mean,
    fun.ymin = mean,
    fun.ymax = mean,
    geom = "crossbar",
    width = 1,
    size = 0.5,
    fatten = 0
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      size = 12,
      colour = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    legend.position = "none",
    strip.text = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  labs(y = "Hue", x = " ") +
  scale_x_discrete(
    labels = c(
      "Col-0",
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
  scale_y_continuous(
    breaks = c(310, 330, 350, 370),
    labels =
      c('310', '330', '350', '10'),
    limits = c(300, 370)
  ) +
  facet_wrap(~ cell.type, ncol = 6) +
  aes(fill = cell) +
  scale_fill_manual(values = barcols) +
  geom_text(
    data = letters.hue.monol,
    aes(label = group),
    family = "Helvetica",
    angle = 90,
    colour = "black",
    hjust = 1,
    size = 3,
    position = position_dodge(width = 0.8)
  )
pdf("hue_monol_rev.pdf", height = 4, width = 10)
b
dev.off()


#### set statistical letters for absorbance ####
phlog.monol.pre %>%
  group_by(cell.type) %>%
  do(data.frame(print.HSD.OD(.)))
letters.OD.monol <-
  read.csv("stats_OD.csv",
           header = FALSE)
colnames(letters.OD.monol) <-
  c("genotype", "mean.OD2", "group", "cell.type")
letters.OD.monol$cell <-
  as.factor(paste(letters.OD.monol$genotype, letters.OD.monol$cell.type))


#### plot absorbance ####
a <- ggplot(phlog.monol.avg, aes(x = genotype, y = mean.OD2)) +
  geom_vline(xintercept = 1,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 3,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 5,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 7,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 9,
             size = 5,
             color = "grey95") +
  geom_vline(xintercept = 11,
             size = 5,
             color = "grey95") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    legend.position = "none",
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
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
    plot.margin = unit(c(0, 0, -0.4, 0), "cm")
  ) +
  labs(y = "Absorbance", x = " ") +
  scale_x_discrete(breaks = c()) +
  scale_y_continuous(
    limits = c(-0.09, 0.55),
    breaks = c(0, 0.2, 0.4),
    labels = c(' 0.0', ' 0.2', ' 0.4')
  ) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.85),
           width = 0.75) +
  facet_wrap(~ cell.type, ncol = 6) +
  aes(fill = cell) +
  scale_fill_manual(values = barcols) +
  geom_jitter(
    data = phlog.monol.pre,
    aes(x = genotype, y = mean.OD1),
    width = 0.1,
    alpha = 0.5
  ) +
  geom_text(
    data = letters.OD.monol,
    aes(label = group),
    family = "Helvetica",
    angle = 90,
    colour = "black",
    hjust = 1,
    size = 3,
    position = position_dodge(width = 0.85)
  )
pdf("OD_monol_rev.pdf", height = 3.7, width = 10)
a
dev.off()

#### import images ####
# WTstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/WT_stained.png"
#     )
#   )
# cl1stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/4cl1_stained.png"
#     )
#   )
# cl2stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/4cl2_stained.png"
#     )
#   )
# cl1x2stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/4cl1x2_stained.png"
#     )
#   )
# ccoaomtstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/ccoaomt_stained.png"
#     )
#   )
# fah1stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/fah1_stained.png"
#     )
#   )
# omt1stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/omt1_stained.png"
#     )
#   )
# ccr1stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/ccr1_stained.png"
#     )
#   )
# cad4stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/cad4_stained.png"
#     )
#   )
# cad5stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/cad5_stained.png"
#     )
#   )
# cad4x5stained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/cad4x5_stained.png"
#     )
#   )
# 
# WTunstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/WT_unstained.png"
#     )
#   )
# cl1unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/4cl1_unstained.png"
#     )
#   )
# cl2unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/4cl2_unstained.png"
#     )
#   )
# cl1x2unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/4cl1x2_unstained.png"
#     )
#   )
# ccoaomtunstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/ccoaomt_unstained.png"
#     )
#   )
# fah1unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/fah1_unstained.png"
#     )
#   )
# omt1unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/omt1_unstained.png"
#     )
#   )
# ccr1unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/ccr1_unstained.png"
#     )
#   )
# cad4unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/cad4_unstained.png"
#     )
#   )
# cad5unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/cad5_unstained.png"
#     )
#   )
# cad4x5unstained <-
#   rasterGrob(
#     readPNG(
#       "~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/cad4x5_unstained.png"
#     )
#   )

#### arrange stained images ####
stained_pre <- plot_grid(
  cl1stained,
  cl2stained,
  cl1x2stained,
  ccoaomtstained,
  fah1stained,
  omt1stained,
  ccr1stained,
  cad4stained,
  cad5stained,
  cad4x5stained,
  labels = c(
    "4cl1",
    "4cl2",
    "4cl1x4cl2",
    "ccoaomt1",
    "fah1",
    "omt1",
    "ccr1",
    "cad4",
    "cad5",
    "cad4xcad5"
  ),
  label_fontfamily = "Helvetica",
  label_fontface = 3,
  scale = 0.98,
  ncol = 5,
  nrow = 2,
  hjust = 0,
  vjust = 1,
  label_x = 0.02,
  label_y = 0.98,
  label_size = 12
)

stained <- plot_grid(
  WTstained,
  stained_pre,
  labels = c("Col-0", ""),
  label_fontfamily = "Helvetica",
  label_fontface = 1,
  scale = 0.99,
  ncol = 2,
  nrow = 1,
  hjust = 0.5,
  vjust = 1,
  label_x = 0.5,
  label_y = 0.98,
  label_size = 12,
  rel_widths = c(2.86, 7.14)
)

#### arrange unstained images ####
unstained_pre <- plot_grid(
  cl1unstained,
  cl2unstained,
  cl1x2unstained,
  ccoaomtunstained,
  fah1unstained,
  omt1unstained,
  ccr1unstained,
  cad4unstained,
  cad5unstained,
  cad4x5unstained,
  labels = c(
    "4cl1",
    "4cl2",
    "4cl1x4cl2",
    "ccoaomt1",
    "fah1",
    "omt1",
    "ccr1",
    "cad4",
    "cad5",
    "cad4xcad5"
  ),
  label_fontfamily = "Helvetica",
  label_fontface = 3,
  scale = 0.98,
  ncol = 5,
  nrow = 2,
  hjust = 0,
  vjust = 1,
  label_x = 0.02,
  label_y = 0.98,
  label_size = 12
)

unstained <- plot_grid(
  WTunstained,
  unstained_pre,
  labels = c("Col-0", ""),
  label_fontfamily = "Helvetica",
  label_fontface = 1,
  scale = 0.99,
  ncol = 2,
  nrow = 1,
  hjust = 0.5,
  vjust = 1,
  label_x = 0.5,
  label_y = 0.99,
  label_size = 12,
  rel_widths = c(2.86, 7.14)
)


#### plot grid ####
pdf("genotypes_grid.pdf", height = 11.75, width = 10)
plot_grid(
  unstained,
  stained,
  a,
  b,
  labels = c('(a)', '(b)', '(c)', '(d)'),
  ncol = 1,
  nrow = 4,
  hjust = 0,
  vjust = 1,
  label_x = 0.001,
  label_y = 0.99,
  label_size = 18,
  label_fontfamily = "Helvetica",
  rel_heights = c(1.2, 1.2, 1, 1.21)
)
dev.off()

#### reduce figure size via ghostscript ####
system(
  "gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r600 -sOutputFile=geno_grid.pdf genotypes_grid.pdf"
)