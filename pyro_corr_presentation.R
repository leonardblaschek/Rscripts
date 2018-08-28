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
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf")
showtext_auto()

py <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/pyrolysis_2018.csv")
py <-
  subset(py, genotype != "pal1" &
           genotype != "pal2" & genotype != "WT_cad")
py$Coniferaldehyde <- py$X22
py$Vanillin <- py$X18
py$Sinapaldehyde <- py$X32
py$Syringaldehyde <- py$X29
py$Aldehydes <-
  ifelse(
    py$variable == "mean",
    py$Coniferaldehyde + py$Sinapaldehyde + py$Vanillin + py$Syringaldehyde,
    sqrt(
      py$Coniferaldehyde ^ 2 + py$Sinapaldehyde ^ 2 + py$Vanillin ^ 2 + py$Syringaldehyde ^ 2
    )
  )
py$Lignin <-
  ifelse(py$variable == "mean", rowSums(py[, 3:28]), sqrt(rowSums(py[, 3:28] ^
                                                                    2)))

py.melt <-
  melt(
    subset(py, select = c(1, 2, 29:34)),
    id = c("genotype", "variable"),
    variable.name = "residue"
  )
py.melt$sd <- py.melt$value
py.melt <-
  merge(
    subset(py.melt, variable == "mean", select = c(1, 3, 4)),
    subset(py.melt, variable == "sd", select = c(1, 3, 5))
  )

py.mean.OD <-
  dcast(subset(phlog.monol.avg, select = c(1, 2, 5)), genotype ~ cell.type)
py.mean.OD[, 7] <-
  0.5 * py.mean.OD[, 2] + 0.21 * py.mean.OD[, 4] + 0.21 * py.mean.OD[, 6] + py.mean.OD[, 3] *
  0.03 + py.mean.OD[, 5] * 0.03
colnames(py.mean.OD)[7] <- "Total"
py.mean.OD <-
  melt(py.mean.OD,
       c("genotype"),
       variable.name = "cell.type",
       value.name = "mean.OD2")
py.mean.SD <-
  dcast(subset(phlog.monol.avg, select = c(1, 2, 6)), genotype ~ cell.type)
py.mean.SD[, 7] <-
  sqrt(
    0.5 * py.mean.SD[, 2] ^ 2 + 0.21 * py.mean.SD[, 4] ^ 2 + 0.21 * py.mean.SD[, 6] ^
      2 + 0.03 * py.mean.SD[, 3] ^ 2 + 0.03 * py.mean.SD[, 5] ^ 2
  )
colnames(py.mean.SD)[7] <- "Total"
py.mean.SD <-
  melt(py.mean.SD,
       c("genotype"),
       variable.name = "cell.type",
       value.name = "SD.OD2")
py.phlog <- merge(py.mean.OD, py.mean.SD)

py.corr <- merge(py.melt, py.phlog)

py.corr$WT <- ifelse(py.corr$genotype == "col-0", "WT", "mutant")

py.corr$residue <- factor(py.corr$residue, levels = c("Coniferaldehyde", "Vanillin", "Sinapaldehyde", "Syringaldehyde", "Aldehydes", "Lignin"))

# plot aldehyde content by genotype
p <-
  ggplot(py.melt,
         aes(
           x = reorder(genotype, -value),
           y = value,
           ymin = value - sd,
           ymax = value + sd
         )) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    colour = "#04253a",
    fill = "grey95"
  ) +
  geom_errorbar(width = 0.2) +
  theme_minimal() +
  facet_wrap(~ residue, ncol = 1, scale = "free_y")

pdf("aldehyde_content.pdf")
p
dev.off()

# calculate linear regressions for log on both axes
file.remove("r_squared.csv")
file.remove("corr.csv")
lin.reg <- function(x) {
  corr.mat <- as.matrix(x[, c(3, 6)])
  corr.mat <- log(corr.mat)
  corr <- data.frame(rcorr(corr.mat)$r)
  corr$p <- data.frame(rcorr(corr.mat)$P)
  corr["cell.type"] <- unique(as.character(x$cell.type))
  corr["residue"] <- unique(as.character(x$residue))
  reg <- lm(log(value) ~ log(mean.OD2), data = x)
  mod <- summary(reg)["adj.r.squared"]
  mod["cell.type"] <- unique(as.character(x$cell.type))
  mod["residue"] <- unique(as.character(x$residue))
  mod["mean.OD2"] <- 2
  mod["value"] <- 0.00001
  mod["SD.OD2"] <- 0
  mod["sd"] <- 0
  write.table(
    corr,
    file = "corr.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
  write.table(
    mod,
    file = "r_squared.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
}

subset(py.corr, cell.type == "Total" & residue != "Aldehydes") %>%
  group_by(cell.type, residue) %>%
  do(data.frame(lin.reg(.)))
adj.r.sq <- read.csv("r_squared.csv", header = FALSE)
adj.r.sq <- adj.r.sq[, 2:8]
colnames(adj.r.sq) <-
  c("r", "cell.type", "residue", "mean.OD2", "value", "SD.OD2", "sd")
adj.r.sq$r <- ifelse(adj.r.sq$r < 0, 0, adj.r.sq$r)

corr <- read.csv("file:///home/leonard/R/Output/wiesner/corr.csv", header = FALSE)
corr <- subset(corr, select = c(1, 2, 4, 6, 7))
colnames(corr) <- c("variable", "r", "p", "cell.type", "residue")
corr <- subset(corr, variable == "mean.OD2", select = c(2:5))
corr["mean.OD2"] <- 1
corr["value"] <- 0.00001
corr["SD.OD2"] <- 0
corr["sd"] <- 0
corr["r"] <-
  ifelse(corr$p < 0.05, paste("r = ", round(corr$r, 3)), "")


# plot log regressions of pyrolysis and wiesner data by cell type and residue
p <-
  ggplot(subset(py.corr, cell.type == "Total" & residue != "Aldehydes"),
         aes(
           x = log(value),
           y = log(mean.OD2),
           ymin = log(mean.OD2) - (0.434 * (SD.OD2 / mean.OD2)),
           ymax = log(mean.OD2) + (0.434 * (SD.OD2 / mean.OD2)),
           xmin = log(value) - (0.434 * (sd / value)),
           xmax = log(value) + (0.434 * (sd / value))
         )) +
  geom_rect(
    data = subset(adj.r.sq),
    aes(fill = r),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.75
  ) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_point(
    aes(colour = WT),
    shape = 16,
    stroke = 0.1,
    size = 2,
    alpha = 0.5
  ) +
  geom_errorbar(aes(colour = WT),
                size = 0.1,
                width = 0.2,
                alpha = 0.8) +
  geom_errorbarh(aes(colour = WT),
                 size = 0.1,
                 height = 0.1,
                 alpha = 0.8) +
  scale_colour_manual(values = c("#04253a", "#04253a")) +
  geom_smooth(
    method = lm,
    colour = "#ffcc3d",
    linetype = 1,
    size = 0.5,
    se = FALSE
  ) +
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
  theme_few() +
  theme(
    text = element_text(size = 12, family = 'Helvetica'),
    strip.text = element_text(
      hjust = 0,
      vjust = 0.5,
      face = "italic",
      colour = "#04253a"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "#04253a", size = 0.25),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 14, colour = "#04253a"),
    axis.title.y = element_text(size = 14, colour = "#04253a"),
    axis.text = element_blank()
  ) +
  labs(
    x = expression(paste("Log"[italic("e")], "(content * biomass" ^ -1, ")")), 
    y = expression(paste("Log"[italic("e")], "(Wiesner Stain)"))
    ) +
  scale_x_continuous(breaks = c(-9, -6, -3)) +
  scale_y_continuous(breaks = c(0, -2, -4), limits = c(-5, 1)) +
  facet_wrap(~ residue, ncol = 5) +
  geom_text(
    data = adj.r.sq,
    aes(label = paste("R² = ", round(r, 3))),
    family = "Helvetica",
    colour = "#04253a",
    size = 3,
    hjust = 0
  ) +
  geom_text(
    data = corr,
    aes(label = r),
    family = "Helvetica",
    colour = "#04253a",
    size = 3,
    hjust = 0
  )

pdf("pyro_corr_both_log_pres.pdf", height = 2, width = 7)
p
dev.off()