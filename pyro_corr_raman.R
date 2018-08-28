library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(colorspace)
library(dplyr)
library(Hmisc)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

raman <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol paper/18-06-25_RAMAN/raman.csv")
raman.avg <- ddply(raman, c("band", "genotype", "cell.type"), summarise,
                   mean = mean(rel.value),
                   sd = sd(rel.value))

raman.mean.OD <-
  dcast(subset(phlog.monol.avg, select = c(1, 2, 5)), genotype ~ cell.type)
raman.mean.OD <-
  melt(raman.mean.OD,
       c("genotype"),
       variable.name = "cell.type",
       value.name = "mean.OD2")
raman.mean.SD <-
  dcast(subset(phlog.monol.avg, select = c(1, 2, 6)), genotype ~ cell.type)
raman.mean.SD <-
  melt(raman.mean.SD,
       c("genotype"),
       variable.name = "cell.type",
       value.name = "SD.OD2")
raman.phlog <- merge(raman.mean.OD, raman.mean.SD)

raman.corr <- merge(raman.avg, subset(raman.phlog, (cell.type == "MX" | cell.type == "XF" | cell.type == "IF") & (genotype == "4cl1x2" | genotype == "col-0")))
# calculate linear regressions for log on both axes
file.remove("r_squared_raman.csv")
file.remove("corr_raman.csv")
lin.reg <- function(x) {
  corr.mat <- as.matrix(x[, c(4, 6)])
  corr <- data.frame(rcorr(corr.mat)$r)
  corr$p <- data.frame(rcorr(corr.mat)$P)
  corr["band"] <- unique(as.character(x$band))
  reg <- lm(mean ~ mean.OD2, data = x)
  mod <- summary(reg)["adj.r.squared"]
  mod["band"] <- unique(as.character(x$band))
  mod["mean.OD2"] <- 0.5
  mod["mean"] <- 0.00001
  mod["SD.OD2"] <- 0
  mod["sd"] <- 0
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
  group_by(band) %>%
  do(data.frame(lin.reg(.)))
adj.r.sq <- read.csv("r_squared_raman.csv", header = FALSE)
adj.r.sq <- adj.r.sq[, 2:7]
colnames(adj.r.sq) <-
  c("r", "band", "mean.OD2", "mean", "SD.OD2", "sd")
adj.r.sq$r <- ifelse(adj.r.sq$r < 0, 0, adj.r.sq$r)

corr <- read.csv("file:///home/leonard/R/Output/wiesner/corr_raman.csv", header = FALSE)
corr <- subset(corr, select = c(1, 2, 4, 6))
colnames(corr) <- c("variable", "r", "p", "band")
corr <- subset(corr, variable == "mean.OD2", select = c(2:4))
corr["mean.OD2"] <- 0.45
corr["mean"] <- 0.00001
corr["SD.OD2"] <- 0
corr["sd"] <- 0
corr["r"] <- paste("r = ", round(corr$r, 3), "| p = ", round(corr$p, 3))


# plot log regressions of ramanrolysis and wiesner data by cell type and residue
p <-
  ggplot(raman.corr,
         aes(
           x = mean,
           y = mean.OD2,
           ymin = mean.OD2 - SD.OD2,
           ymax = mean.OD2 + SD.OD2,
           xmin = mean - sd,
           xmax = mean + sd
         )) +
  geom_rect(
    data = subset(adj.r.sq),
    aes(fill = r),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.5
  ) +
  geom_smooth(
    method = lm,
    colour = "black",
    linetype = 2,
    size = 0.5,
    se = FALSE
  ) +
  scale_fill_distiller(palette = "Greens", direction = 1, guide = FALSE) +
  geom_errorbar(size = 0.25,
                width = 0.1,
                alpha = 0.8) +
  geom_errorbarh(size = 0.25,
                 height = 0.01,
                 alpha = 0.8) +
  geom_point(aes(colour = cell.type, shape = genotype),
    stroke = 0.1,
    size = 2,
    alpha = 1
  ) +
  scale_color_brewer(palette = "Dark2") +

  # geom_ribbon(
  #   stat = "smooth",
  #   method = "lm",
  #   fill = "blue",
  #   linetype = 2,
  #   colour = NA,
  #   size = 0.2,
  #   alpha = 0.1,
  #   level = 0.95
  # ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = 'Helvetica'),
    strip.text = element_text(
      hjust = 0,
      vjust = 0.5,
      face = "italic",
      colour = "black"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.ticks.x = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      colour = "black"
    )
  ) +
  labs(
    x = "Raman intensity",
    y = "Wiesner absorbance"
    ) +
  # scale_x_continuous(breaks = c(-9, -6, -3)) +
  # scale_y_continuous(breaks = c(0, -2, -4)) +
  facet_wrap(~ band, ncol = 3) +
  geom_text(
    data = adj.r.sq,
    aes(label = paste("R² = ", round(r, 3))),
    family = "Helvetica",
    colour = "black",
    size = 3,
    hjust = 0
  ) +
  geom_text(
    data = corr,
    aes(label = r),
    family = "Helvetica",
    colour = "black",
    size = 3,
    hjust = 0
  )

pdf("ramanro_corr_both_log.pdf", height = 6)
p
dev.off()