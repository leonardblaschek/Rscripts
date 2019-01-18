library(sysfonts)
library(showtext)
library(ggplot2)
library(gplots)
library(ggthemes)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(colorspace)
library(dplyr)
library(Hmisc)
library(corrplot)
library(cowplot)

#### import Helvetica Neue ####
font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

#### import measurements for Wiesner absorbance in sections ####
phlog.monol <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
           skip = 2)

#### calculate pixel values from OD ####
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


#### calculate the correct hue on the 360 point circular scale ####
phlog.monol$hue <- ((phlog.monol$h.stained + 128) / 255 * 360)

phlog.monol$replicate <-
  as.factor(as.character(phlog.monol$replicate))


#### calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem ####
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

#### average per replicate (for boxplots) ####
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


#### average per genotype (for barplots) ####
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

#### import Raman measurements ####
rspec <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/18-06-25_RAMAN/raman_spectra.csv"
  )

raman.phlog <-
  dcast(subset(phlog.monol.avg, select = c(1, 2, 5)), genotype ~ cell.type)
raman.phlog <-
  melt(raman.phlog,
       c("genotype"),
       variable.name = "cell.type",
       value.name = "absorbance")

raman.corr <-
  merge(rspec, subset(
    raman.phlog,
    (cell.type == "MX" |
       cell.type == "XF" |
       cell.type == "IF") &
      (genotype == "4cl1x2" | genotype == "col-0")
  ))

#### calculate linear regressions ####
file.remove("r_squared_raman.csv")
file.remove("corr_raman.csv")
lin.reg <- function(x) {
  corr.mat <- as.matrix(x[, c(4, 5)])
  corr <- data.frame(rcorr(corr.mat)$r)
  corr$p <- data.frame(rcorr(corr.mat)$P)
  corr["wave.number"] <- unique(as.character(x$wave.number))
  reg <- lm(mean ~ absorbance, data = x)
  mod <- summary(reg)["adj.r.squared"]
  mod["wave.number"] <- unique(as.character(x$wave.number))
  write.table(
    corr,
    file = "corr_raman.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
  write.table(
    mod,
    file = "r_squared_raman.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
}

raman.corr %>%
  group_by(wave.number) %>%
  do(data.frame(lin.reg(.)))

adj.r.sq <- read.csv("r_squared_raman.csv", header = FALSE)
adj.r.sq <- adj.r.sq[, 2:3]
adj.r.sq <- adj.r.sq[!duplicated(adj.r.sq),]
colnames(adj.r.sq) <-
  c("r2", "wave.number")

corr <-
  read.csv("corr_raman.csv",
           header = FALSE)
corr <- subset(corr, select = c(1, 2, 4, 6))
corr <- corr[!duplicated(corr),]
colnames(corr) <- c("variable", "r", "p", "wave.number")
corr <- subset(corr, variable == "absorbance", select = c(2:4))
corr$y <- 1

corr <- merge(corr, adj.r.sq)

#### plot figure (needs rearrangement in inkscape) ####
pdf("raman_corr.pdf")
corr.heat <-
  ggplot(data = corr, aes(x = wave.number, y = y, fill = r)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c(
      "#2166ac",
      "#67a9cf",
      "#d1e5f0",
      "white",
      "#fddbc7",
      "#ef8a62",
      "#b2182b"
    ),
    limits = c(-1, 1)
  ) +
  theme_void() +
  scale_x_reverse() +
  theme(
    legend.position = c(0.5, 0.5),
    legend.key.height = unit(1.5, "mm"),
    legend.key.width = unit(3.5, "mm"),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(family = "Helvetica", size = 6),
    plot.margin = unit(c(-1.1, 0, 1.12, 0), "cm")
  )

p.heat <-
  ggplot(data = corr, aes(x = wave.number, y = y, fill = p)) +
  geom_tile() +
  scale_fill_gradientn(
    trans = "log",
    limits = c(0.01, 1),
    breaks = c(0.01, 0.05, 0.1, 0.5),
    colours = c("#FDE725FF", "#73D055FF", "#238A8DFF", "#440154FF")
  ) +
  theme_void() +
  scale_x_reverse() +
  theme(
    legend.position = c(0.5, 0.5),
    legend.key.height = unit(1.5, "mm"),
    legend.key.width = unit(3.5, "mm"),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(family = "Helvetica", size = 6),
    plot.margin = unit(c(-1.1, 0, 1.12, 0), "cm")
  )

rspec.col.MX <-
  ggplot(data = subset(rspec, cell.type == "MX" &
                         genotype == "col-0"),
         aes(x = wave.number, y = mean)) +
  geom_vline(xintercept = 1596.4,
             colour = "grey75",
             size = 0.5) +
  geom_vline(
    xintercept = 1620,
    colour = "grey75",
    size = 0.5,
    linetype = 2
  ) +
  geom_vline(
    xintercept = 1140,
    colour = "grey75",
    size = 0.5,
    linetype = 2
  ) +
  geom_line(group = 1, size = 0.5) +
  theme_few() +
  expand_limits(y = -150) +
  scale_x_reverse(breaks = c(750, 1000, 1250, 1500, 1750)) +
  labs(x = expression(paste('Wave number [cm' ^ {
    '-1'
  }, ']')),
  y = "Intensity") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    axis.title.y = element_text(size = 14, margin = margin(
      t = 0,
      r = -1,
      b = 0,
      l = 0
    )),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = "black", size = 14),
    axis.ticks.x = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    plot.margin = unit(c(0.75, 0.75, -1.25, 0.75), "cm")
  )
dev.off()

pdf("rspec_grid.pdf", height = 3, width = 8)
plot_grid(
  rspec.col.MX,
  corr.heat,
  p.heat,
  ncol = 1,
  rel_heights = c(1, 0.1, 0.1),
  rel_widths = c(1, 1, 1),
  align = "v"
)
dev.off()