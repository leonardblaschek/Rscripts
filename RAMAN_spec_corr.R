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

font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

rspec <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/18-06-25_RAMAN/raman_spectra.csv"
  )
rspec_colMX <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/18-06-25_RAMAN/WT_MX_1.csv")

rspec.norm <- subset(rspec, wave.number == 1599.8)
rspec.norm[, 5] <- subset(rspec, wave.number == 1095.3, select = 4)
colnames(rspec.norm)[4] <- "lignin"
colnames(rspec.norm)[5] <- "cellulose"

rspec <-
  merge(rspec,
        rspec.norm[, 2:5],
        by = c("genotype", "cell.type"),
        all = TRUE)

# CALCULATE RATIOS
# rspec$mean <- rspec$mean / rspec$lignin

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
       cell.type == "IF") & (genotype == "4cl1x2" | genotype == "col-0")
  ))

# calculate linear regressions for log on both axes
file.remove("r_squared_raman.csv")
file.remove("corr_raman.csv")
lin.reg <- function(x) {
  corr.mat <- as.matrix(x[, c(4, 7)])
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
adj.r.sq <- adj.r.sq[!duplicated(adj.r.sq), ]
colnames(adj.r.sq) <-
  c("r2", "wave.number")

corr <-
  read.csv("file:///home/leonard/R/Output/wiesner/corr_raman.csv",
           header = FALSE)
corr <- subset(corr, select = c(1, 2, 4, 6))
corr <- corr[!duplicated(corr), ]
colnames(corr) <- c("variable", "r", "p", "wave.number")
corr <- subset(corr, variable == "absorbance", select = c(2:4))
corr$y <- 1

corr <- merge(corr, adj.r.sq)


pdf("raman_corr.pdf")
corr.heat <-
  ggplot(data = corr, aes(x = wave.number, y = y, fill = r)) +
  geom_tile() +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  # expand_limits(fill=1) +
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
corr.heat

p.heat <-
  ggplot(data = corr, aes(x = wave.number, y = y, fill = p)) +
  geom_tile() +
  # scale_fill_viridis_c(trans = "log", direction = -1, limits = c(0,1)) +
  scale_fill_gradientn(
    trans = "log",
    limits = c(0.01, 1),
    breaks = c(0.01, 0.05, 0.1, 0.5),
    colours = c("#FDE725FF", "#73D055FF","#238A8DFF", "#440154FF")
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
p.heat

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
  # scale_y_continuous(limits = c(0, 1000)) +
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
rspec.col.MX
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

# rspec_mod <- rasterGrob(readPNG("/home/leonard/Documents/LaTex/figures_edouard/RAMAN_corr.pdf"))
# pdf("rspec_heat.pdf", height = 10, width = 6)
# plot_grid(fire,
#           rspec_mod,
#           labels = c("", "M"),
#           label_fontfamily = "Helvetica",
#           label_colour = "white",
#           scale = 0.98,
#           ncol = 1,
#           nrow = 2,
#           hjust = 0,
#           vjust = 1,
#           label_x = 0.02,
#           label_y = 0.98,
#           rel_heights = c(0.8, 0.4))
# dev.off()
# pdf("rspec_grid_complete.pdf", height = 6, width = 5)
# plot_grid(
#   raw,
#   cellu,
#   lig,
#   labels = c("Absolute", "Norm. to cellulose", "Norm. to lignin"),
#   label_fontfamily = "Helvetica",
#   label_fontface = "plain",
#   ncol = 1,
#   align = "v",
#   hjust = 0,
#   label_x = 0.01
# )
# dev.off()
