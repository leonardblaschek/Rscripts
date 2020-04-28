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

#### import Helvetica Neue ####
font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
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

#### import pyrolysis data, and summarise peaks ####
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

write_csv(py, "pyrolysis_At.csv")

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

#### merge absorbance and pyrolysis data ####
py.mean.OD <-
  dcast(subset(phlog.monol.avg, genotype != "ccr1xfah1", select = c(1, 2, 5)),
        genotype ~ cell.type)
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

py.corr$residue <-
  factor(
    py.corr$residue,
    levels = c(
      "Coniferaldehyde",
      "Vanillin",
      "Sinapaldehyde",
      "Syringaldehyde",
      "Aldehydes",
      "Lignin"
    )
  )

#### calculate log-log linear regressions ####
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

#### calculate pearson correlations ####
py.corr %>%
  group_by(cell.type, residue) %>%
  do(data.frame(lin.reg(.)))
adj.r.sq <- read.csv("r_squared.csv", header = FALSE)
adj.r.sq <- adj.r.sq[, 2:8]
colnames(adj.r.sq) <-
  c("r", "cell.type", "residue", "mean.OD2", "value", "SD.OD2", "sd")
adj.r.sq$r <- ifelse(adj.r.sq$r < 0, 0, adj.r.sq$r)

corr <- read.csv("corr.csv", header = FALSE)
corr <- subset(corr, select = c(1, 2, 4, 6, 7))
colnames(corr) <- c("variable", "r", "p", "cell.type", "residue")
corr <- subset(corr, variable == "mean.OD2", select = c(2:5))
corr["mean.OD2"] <- 1
corr["value"] <- 0.00001
corr["SD.OD2"] <- 0
corr["sd"] <- 0
corr["r"] <-
  ifelse(corr$p < 0.05, paste("r = ", round(corr$r, 3)), "")


#### plot log regressions of pyrolysis and wiesner data by cell type and residue ####
# logarithmic error bars according to https://faculty.washington.edu/stuve/log_error.pdf
p <-
  ggplot(py.corr,
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
    alpha = 0.5
  ) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  geom_point(
    aes(colour = WT),
    shape = 16,
    stroke = 0.1,
    size = 1,
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
  scale_colour_manual(values = c("#000000", "#e20000")) +
  geom_smooth(
    method = lm,
    colour = "blue",
    linetype = 1,
    size = 0.25,
    se = FALSE
  ) +
  theme_few() +
  theme(
    text = element_text(size = 12, family = 'Helvetica'),
    strip.text.x = element_text(
      size = 6,
      hjust = 0,
      vjust = 0,
      face = "italic",
      colour = "black"
    ),
    strip.text.y = element_text(
      size = 6,
      hjust = 0.5,
      vjust = 0,
      face = "italic",
      colour = "black"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title = element_text(size = 6, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
  ) +
  labs(x = expression(paste("Log"[italic("e")], "(content * biomass" ^ -1, ")")),
       y = expression(paste("Log"[italic("e")], "(wiesner stain absorbance)"))) +
  scale_x_continuous(breaks = c(-9, -6, -3)) +
  scale_y_continuous(breaks = c(0, -2, -4), limits = c(-5, 1)) +
  facet_grid(residue ~ cell.type) +
  geom_text(
    data = adj.r.sq,
    aes(label = paste("R² = ", round(r, 3))),
    family = "Helvetica",
    colour = "black",
    size = 6 / (14/5),
    hjust = 0
  ) +
  geom_text(
    data = corr,
    aes(label = r),
    family = "Helvetica",
    colour = "black",
    size = 6 / (14/5),
    hjust = 0
  )

pdf("pyro_corr_both_log.pdf", width = 5, height = 5)
p
dev.off()

#### papperstidning figure ####

theme_leo <- function(base_size = 8,
                      base_family = "Helvetica") {
  theme_minimal(
    base_size = base_size,
    base_family = base_family
  ) %+replace%
    theme(
      strip.text = element_text(hjust = 0, face = "italic"),
      axis.ticks = element_line(
        size = 0.2,
        lineend = "square",
        color = "black"
      ),
      axis.text.x = element_text(
        size = 8,
        colour = "black", # flipped coords
        margin = margin(1, 1, 1, 1)
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 8,
        angle = 0,
        vjust = 0.5,
        hjust = 1,
        margin = margin(1, 1, 1, 1)
      ),
      axis.title = element_text(
        colour = "black",
        size = 8
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.2),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "none",
      legend.text = element_text(size = 8),
      legend.key.height = unit(4, "mm"),
      legend.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.margin = unit(c(2, 2, 2, 2), "mm"),
      complete = TRUE
    )
}

py_papper <- subset(py.corr, residue == "Coniferaldehyde" & cell.type == "Total")

papper_plot <- ggplot(py_papper,
       aes(
         x = log(value),
         y = log(mean.OD2),
         ymin = log(mean.OD2) - (0.434 * (SD.OD2 / mean.OD2)),
         ymax = log(mean.OD2) + (0.434 * (SD.OD2 / mean.OD2)),
         xmin = log(value) - (0.434 * (sd / value)),
         xmax = log(value) + (0.434 * (sd / value))
       )) +
  # geom_rect(
  #   data = subset(adj.r.sq),
  #   # aes(fill = r),
  #   xmin = -Inf,
  #   xmax = Inf,
  #   ymin = -Inf,
  #   ymax = Inf,
  #   alpha = 0.5
  # ) +
  # scale_fill_distiller(palette = "Greens", direction = 1) +
  geom_point(
    # aes(colour = WT),
    shape = 16,
    size = 2,
    alpha = 0.5
  ) +
  geom_errorbar(
    # aes(colour = WT),
                size = 0.2,
                width = 0.05,
                alpha = 0.8) +
  geom_errorbarh(
    # aes(colour = WT),
                 size = 0.2,
                 height = 0.025,
                 alpha = 0.8) +
  # scale_colour_manual(values = c("#000000", "#e20000")) +
  geom_smooth(
    method = lm,
    colour = "black",
    linetype = 2,
    size = 0.25,
    se = FALSE
  ) +
  theme_leo() +
  # theme(
  #   text = element_text(size = 12, family = 'Helvetica'),
  #   strip.text.x = element_text(
  #     size = 6,
  #     hjust = 0,
  #     vjust = 0,
  #     face = "italic",
  #     colour = "black"
  #   ),
  #   strip.text.y = element_text(
  #     size = 6,
  #     hjust = 0.5,
  #     vjust = 0,
  #     face = "italic",
  #     colour = "black"
  #   ),
  #   panel.grid.minor = element_blank(),
  #   panel.grid.major.x = element_blank(),
  #   legend.position = "none",
  #   panel.border = element_rect(fill = NA, color = "black", size = 0.25),
  #   axis.ticks = element_line(
  #     size = 0.25,
  #     lineend = "square",
  #     color = "black"
  #   ),
  #   axis.title = element_text(size = 6, colour = "black"),
  #   axis.text = element_text(size = 6, colour = "black"),
  # ) +
  labs(x = expression(paste("Log"[italic("e")], "(coniferaldehyd)")),
       y = expression(paste("Log"[italic("e")], "(wiesner test)")))
  # scale_x_continuous(breaks = c(-9, -6, -3)) +
  # scale_y_continuous(breaks = c(0, -2, -4), limits = c(-5, 1)) +
  # facet_grid(residue ~ cell.type) 
  # geom_text(
  #   data = subset(adj.r.sq, residue == "Coniferaldehyde" & cell.type == "Total"),
  #   aes(label = paste("R² = ", round(r, 3))),
  #   y = -1,
  #   family = "Helvetica",
  #   colour = "black",
  #   size = 6 / (14/5),
  #   hjust = 0
  # ) +
  # geom_text(
  #   data = subset(corr, residue == "Coniferaldehyde" & cell.type == "Total"),
  #   aes(label = r),
  #   y = -1.2,
  #   family = "Helvetica",
  #   colour = "black",
  #   size = 6 / (14/5),
  #   hjust = 0
  # )
pdf("papper_plot.pdf", height = 1.5, width = 2)
papper_plot
dev.off()

#### supplemental regressions ####

py <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/pyrolysis_2018.csv")
py <-
  subset(py, genotype != "pal1" &
           genotype != "pal2" & genotype != "WT_cad")
py$H <-
  ifelse(py$variable == "mean", rowSums(py[, 7:9]), sqrt(rowSums(py[, 7:9]^
                                                                   2)))
py$G <-
  ifelse(py$variable == "mean", rowSums(py[, 7:19]), sqrt(rowSums(py[, 7:19]^
                                                                    2)))
py$S <-
  ifelse(py$variable == "mean", rowSums(py[, 21:28]), sqrt(rowSums(py[, 21:28]^
                                                                     2)))
py.melt <-
  melt(
    subset(py, select = c(1, 2, 30:32)),
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
  dcast(subset(phlog.monol.avg, genotype != "ccr1xfah1", select = c(1, 2, 5)),
        genotype ~ cell.type)
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


supp_reg <- py.corr %>%
  filter(cell.type == "Total") %>%
  pivot_wider(id_cols = c(genotype, mean.OD2), names_from = residue, values_from = value)

mod_G <- lm(mean.OD2 ~ G, data = supp_reg)
mod_GS <- lm(mean.OD2 ~ G + S, data = supp_reg)
mod_GH <- lm(mean.OD2 ~ G + H, data = supp_reg)
mod_GSH <- lm(mean.OD2 ~ G + S + H, data = supp_reg)

models <- broom::glance(mod_G) %>%
  rbind(broom::glance(mod_GS))%>%
  rbind(broom::glance(mod_GH))%>%
  rbind(broom::glance(mod_GSH)) %>%
  add_column(model = c("G", "G + S", "G + H", "G + S + H"), .before = 1) %>%
  select(model, adj.r.squared, BIC)

model_table <- kable(models,
                      "latex",
                      align = "lll",
                      caption = "Multiple regression models of the Wiesner test dependency on the concentration of different lignin subunits.",
                      col.names = c(
                        "Predictors",
                        "Adjusted R²",
                        "BIC"
                      ),
                      booktabs = TRUE,
                      escape = TRUE,
                      linesep = ""
) 
# %>%
#   add_header_above(c(
#     " " = 1, "Relative AUC\\\\textsuperscript{a}" = 2,
#     "Signal\\\\textsubscript{Raman} * Signal\\\\textsubscript{fluo}\\\\textsuperscript{-1}" = 2
#   ),
#   escape = FALSE
#   ) %>%
#   footnote(alphabet = "Relative to the AUC of monomeric G\\\\textsubscript{CHOH}", escape = FALSE)

model_table %>%
  save_kable("Wiesner_tableS2.pdf")

readr::write_file(model_table, "model_table.txt")
