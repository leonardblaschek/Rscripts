library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(showtext)
library(colorscience)
library(cowplot)
library(tidyverse)
data(illuminants)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

peaks <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-11-15_spectras_platereader/17-11-15_collected_peaks.csv"
  )
peaks <- melt(peaks, id = c("peak", "repeat.", "condition"))
peaks$variable <- gsub("X([0-9]+)", "\\1", peaks$variable)
peaks$variable <- as.numeric(peaks$variable)
peaks$peak <- (paste('Peak', peaks$peak))
peaks$condition <-
  revalue(peaks$condition,
          c('stained' = 'acidified', 'unstained' = 'untreated'))
peaks$peak <-
  ordered(
    peaks$peak,
    levels = c(
      "Peak 1",
      "Peak 2",
      "Peak 3",
      "Peak 4",
      "Peak 5",
      "Peak 6",
      "Peak 7",
      "Peak 8",
      "Peak 9",
      "Peak 10",
      "Peak 11",
      "Peak 12"
    )
  )
peaks.lbl <- peaks[1:12,]
peaks.lbl$variable <- 280
peaks.lbl$value <- 0.3
peaks <- subset(peaks, peaks$peak == "Peak 1" | peaks$peak == "Peak 5" | peaks$peak == "Peak 6" | peaks$peak == "Peak 9" | peaks$peak == "Peak 12")
peaks.lbl <- subset(peaks.lbl, peaks.lbl$peak == "Peak 1" | peaks.lbl$peak == "Peak 5" | peaks.lbl$peak == "Peak 6" | peaks.lbl$peak == "Peak 9" | peaks.lbl$peak == "Peak 12")
# peaks$peak <- as.factor(peaks$peak)
# peaks.avg <- ddply(peaks, c("peak","condition", "variable"), summarise, mean=mean(value))

# # calculate HSV colours for Peak 6/9
spec.peak9 <-
  subset(
    peaks,
    peak == "Peak 9" &
      repeat. == 1 &
      condition == "acidified" &
      variable > 359 & variable < 800,
    select = c(4, 5)
  )
spec.peak9 <- spec.peak9[seq(1, nrow(spec.peak9), 5), ]
spec.peak9$value <-  (10 ^ -spec.peak9$value) * 100
spec.peak9 <- data.matrix(spec.peak9)
XYZ.peak9 <-
  spectra2XYZ(spec.peak9, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.peak9 <- XYZ2RGB(XYZ.peak9, illuminant = 'E')
HSV.peak9 <- RGB2HSV(RGB.peak9)
hue.peak9 <- HSV.peak9[, 1] * 360
names(hue.peak9) <- NULL

peak9.col <- data.frame("Peak 9", "acidified", 1, "p9", 280, 0.225, round(hue.peak9, 0))
colnames(peak9.col) <- c("peak", "condition", "repeat.", "fll", "variable", "value", "lbl")

pdf("spectra_peaks.pdf", width = 6, height = 1.5)
SP <-
  ggplot(subset(peaks, repeat. == 1),
         aes(y = value, x = variable, linetype = condition)) +
  geom_rect(
    data = peak9.col,
    aes(fill = fll),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = c("p9" = hsv(HSV.peak9[, 1], 10*HSV.peak9[, 2], 1)),
    guide = FALSE
  ) +
  # geom_vline(xintercept = 550,
  #            color = "grey35",
  #            linetype = 1) +
  geom_line() +
  theme_minimal() +
  theme(
    text = element_text(size = 15, family = "Helvetica"),
    strip.text = element_blank(),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      vjust = 1,
      hjust = 1,
      color = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing = unit(1.5, "mm"),
    # plot.margin = unit(c(0, 0, 0, 0), "cm"),
    # legend.position = c(0.12, 0.86),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.height = unit(4, "mm")
  ) +
  xlim(270, 650) +
  #     ylim(0, 0.5) +
  scale_y_continuous(breaks = c(0, 0.2), limits = c(0, 0.35), labels = c("0", "0.2")) +
  labs(x = "Wavelength [nm]", y = "Absorbance") +
  facet_wrap(~ peak, ncol = 5) +
  geom_text(
    data = peaks.lbl,
    aes(label = peak),
    family = "Helvetica",
    colour = "black",
    fontface = "italic",
    size = 3,
    hjust = 0) +
  geom_text(
      data = peak9.col,
      aes(label = paste("calc. hue: ", lbl)),
      family = "Helvetica",
      colour = "black",
      size = 3,
      hjust = 0
    ) 
SP 
dev.off()

monomers50 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-11-15_spectras_platereader/17-11-15_monomers_50uM.csv"
  )
monomers50 <-
  reshape2::melt(monomers50, id = c("Monomer", "Acid", "Time", "Concentration"))
monomers50$variable <- gsub("X([0-9]+)", "\\1", monomers50$variable)
monomers50$variable <- as.numeric(monomers50$variable)
monomers50$Monomer <-
  revalue(monomers50$Monomer, c('Conifer alcohol' = 'Coniferyl alcohol'))
monomers50$Monomer <-
  ordered(
    monomers50$Monomer,
    levels = c(
      "Coniferyl alcohol",
      "Coniferaldehyde",
      "Sinapaldehyde",
      "Vanillin",
      "Syringaldehyde"
    )
  )

# calculate HSV colours for G-CHO at 250 mM
spec.mGCHO250 <-
  subset(
    monomers50,
    Monomer == "Coniferaldehyde" &
      Acid == "250 mM HCl" &
      variable > 359 & variable < 800 & Time == "0 h",
    select = c(5, 6)
  )
spec.mGCHO250 <- spec.mGCHO250[seq(1, nrow(spec.mGCHO250), 5), ]
spec.mGCHO250$value <-  (10 ^ -spec.mGCHO250$value) * 100
spec.mGCHO250 <- data.matrix(spec.mGCHO250)
XYZ.mGCHO250 <-
  spectra2XYZ(spec.mGCHO250, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.mGCHO250 <- XYZ2RGB(XYZ.mGCHO250, illuminant = 'E')
HSV.mGCHO250 <- RGB2HSV(RGB.mGCHO250)
hue.mGCHO250 <- HSV.mGCHO250[, 1] * 360
names(hue.mGCHO250) <- NULL

# calculate HSV colours for G-CHO at 3 M
spec.mGCHO3 <-
  subset(
    monomers50,
    Monomer == "Coniferaldehyde" &
      Acid == "3 M HCl" &
      variable > 359 & variable < 800 & Time == "0 h",
    select = c(5, 6)
  )
spec.mGCHO3 <- spec.mGCHO3[seq(1, nrow(spec.mGCHO3), 5), ]
spec.mGCHO3$value <-  (10 ^ -spec.mGCHO3$value) * 100
spec.mGCHO3 <- data.matrix(spec.mGCHO3)
XYZ.mGCHO3 <-
  spectra2XYZ(spec.mGCHO3, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.mGCHO3 <- XYZ2RGB(XYZ.mGCHO3, illuminant = 'E')
HSV.mGCHO3 <- RGB2HSV(RGB.mGCHO3)
hue.mGCHO3 <- HSV.mGCHO3[, 1] * 360
names(hue.mGCHO3) <- NULL

# calculate HSV colours for S-CHO at 250 mM
spec.mSCHO250 <-
  subset(
    monomers50,
    Monomer == "Sinapaldehyde" &
      Acid == "250 mM HCl" &
      variable > 359 & variable < 800 & Time == "0 h",
    select = c(5, 6)
  )
spec.mSCHO250 <- spec.mSCHO250[seq(1, nrow(spec.mSCHO250), 5), ]
spec.mSCHO250$value <-  (10 ^ -spec.mSCHO250$value) * 100
spec.mSCHO250 <- data.matrix(spec.mSCHO250)
XYZ.mSCHO250 <-
  spectra2XYZ(spec.mSCHO250, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.mSCHO250 <- XYZ2RGB(XYZ.mSCHO250, illuminant = 'E')
HSV.mSCHO250 <- RGB2HSV(RGB.mSCHO250)
hue.mSCHO250 <- HSV.mSCHO250[, 1] * 360
names(hue.mSCHO250) <- NULL

# calculate HSV colours for S-CHO at 3 mM
spec.mSCHO3 <-
  subset(
    monomers50,
    Monomer == "Sinapaldehyde" &
      Acid == "3 M HCl" &
      variable > 359 & variable < 800 & Time == "0 h",
    select = c(5, 6)
  )
spec.mSCHO3 <- spec.mSCHO3[seq(1, nrow(spec.mSCHO3), 5), ]
spec.mSCHO3$value <-  (10 ^ -spec.mSCHO3$value) * 100
spec.mSCHO3 <- data.matrix(spec.mSCHO3)
XYZ.mSCHO3 <-
  spectra2XYZ(spec.mSCHO3, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.mSCHO3 <- XYZ2RGB(XYZ.mSCHO3, illuminant = 'E')
HSV.mSCHO3 <- RGB2HSV(RGB.mSCHO3)
hue.mSCHO3 <- HSV.mSCHO3[, 1] * 360
names(hue.mSCHO3) <- NULL

# labels for plot
labels <-
  data.frame("Coniferaldehyde",
             "250 mM HCl",
             "0 h",
             "50 µM",
             280,
             3,
             round(hue.mGCHO250, 0),
             "mGCHO250")
labels[, 1] <-
  factor("Coniferaldehyde",
         levels = c("Coniferaldehyde", "Sinapaldehyde"))
labels[, 2] <-
  factor("250 mM HCl", levels = c("250 mM HCl", "3 M HCl"))
labels[, 8] <-
  factor("mGCHO250",
         levels = c("mGCHO250", "mGCHO3", "mSCHO250", "mSCHO3"))
labels[2,] <-
  c("Coniferaldehyde",
    "3 M HCl",
    "0 h",
    "50 µM",
    280,
    3,
    round(hue.mGCHO3, 0),
    "mGCHO3")
labels[3,] <-
  c("Sinapaldehyde",
    "250 mM HCl",
    "0 h",
    "50 µM",
    280,
    3,
    round(hue.mSCHO250, 0),
    "mSCHO250")
labels[4,] <-
  c("Sinapaldehyde",
    "3 M HCl",
    "0 h",
    "50 µM",
    280,
    3,
    round(hue.mSCHO3, 0),
    "mSCHO3")
colnames(labels) <-
  c("Monomer",
    "Acid",
    "Time",
    "Concentration",
    "variable",
    "value",
    "lbl",
    "fll")
labels$variable <- as.numeric(as.character(labels$variable))
labels$value <- as.numeric(as.character(labels$value))


pdf("spectra_monomers.pdf")
MSL <-
  ggplot(monomers50, aes(y = value, x = variable, group = Time)) +
  geom_rect(
    data = labels,
    aes(fill = fll),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "mGCHO250" = hsv(HSV.mGCHO250[, 1], HSV.mGCHO250[, 2], 1),
      "mGCHO3" = hsv(HSV.mGCHO3[, 1], HSV.mGCHO3[, 2], 1),
      "mSCHO250" = hsv(HSV.mSCHO250[, 1], HSV.mSCHO250[, 2], 1),
      "mSCHO3" = hsv(HSV.mSCHO3[, 1], HSV.mSCHO3[, 2], 1)),
      guide = FALSE
    ) +
    geom_vline(xintercept = 561,
               color = "grey35",
               linetype = 1) +
    geom_vline(xintercept = 550,
               color = "grey35",
               linetype = 1) +
    annotate(
      "text",
      x = 567,
      y = 2.2,
      label = '561',
      size = 3,
      hjust = 0,
      color = "grey35",
      family = 'Helvetica'
    ) +
    annotate(
      "text",
      x = 547,
      y = 2.2,
      label = '550',
      size = 3,
      hjust = 1,
      color = "grey35",
      family = 'Helvetica'
    ) +
    geom_line(aes(linetype = Time)) +
      scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(0, 3.1)) +
      theme_minimal() +
      theme(
        text = element_text(size = 15, family = "Helvetica"),
        strip.text = element_text(hjust = 0, face = "italic"),
        # axis.line.y = element_line(size = 0.75, lineend = "square"),
        axis.ticks = element_line(
          size = 0.25,
          lineend = "square",
          color = "grey35"
        ),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(
          size = 10,
          angle = 0,
          vjust = 1,
          hjust = 0.5
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey35", size = 0.25),
        panel.spacing = unit(1.5, "mm"),
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = c(0.375, 0.975),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.key.height = unit(4, "mm")
      ) +
      xlim(270, 650) +
      #     ylim(0, 1) +
      labs(x = "Wavelength [nm]", y = "Absorbance") +
  geom_text(
    data = labels,
    aes(label = paste("calc. hue: ", lbl)),
    family = "Helvetica",
    colour = "black",
    size = 4,
    hjust = 0
  ) +
      facet_grid(Monomer ~ Acid)
    MSL
    dev.off()
    
    pdf("spectra_grid.pdf", height = 8, width = 8)
    right_col <- plot_grid(
      DSLS,
      DSS,
      labels = c('B', 'C'),
      ncol = 1,
      nrow = 2,
      label_fontfamily = "Helvetica"
    )
    left_col <- plot_grid(
      MSL,
      SP,
      labels = c('A', 'B'),
      ncol = 1,
      nrow = 2,
      rel_heights = c(2.5,1),
      label_fontfamily = "Helvetica"
    )
    plot_grid(
      MSL,
      right_col,
      labels = c('A', ''),
      label_fontfamily = "Helvetica"
    )
    dev.off()
    
rev_spectra <- read_csv("/home/leonard/Documents/Uni/Phloroglucinol/2019-12-20_spectra_platereader/2019-12-20_spectra_stained.csv") %>%
  pivot_longer(cols = -c("sample", "condition"), names_to = "variable", values_to = "value") %>%
  mutate(variable = as.numeric(variable))

plotly::ggplotly(ggplot(rev_spectra, aes(x = variable, y = value, colour = sample)) + geom_line())

spec.Hcho <-
  subset(
    rev_spectra,
    sample == "Hcho" &
      variable > 359 & variable < 800,
    select = c(3, 4)
  )
spec.Hcho <- spec.Hcho[seq(1, nrow(spec.Hcho), 5), ]
spec.Hcho$value <-  (10 ^ -spec.Hcho$value) * 100
spec.Hcho <- data.matrix(spec.Hcho)
XYZ.Hcho <-
  spectra2XYZ(spec.Hcho, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.Hcho <- XYZ2RGB(XYZ.Hcho, illuminant = 'E')
HSV.Hcho <- RGB2HSV(RGB.Hcho)
hue.Hcho <- HSV.Hcho[, 1] * 360
names(hue.Hcho) <- NULL

spec.Gdhp <-
  subset(
    rev_spectra,
    sample == "Gdhp" &
      variable > 359 & variable < 800,
    select = c(3, 4)
  )

spec.Gdhp <- spec.Gdhp[seq(1, nrow(spec.Gdhp), 5), ]
spec.Gdhp$value <-  (10 ^ -spec.Gdhp$value) * 100
spec.Gdhp <- data.matrix(spec.Gdhp)
XYZ.Gdhp <-
  spectra2XYZ(spec.Gdhp, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.Gdhp <- XYZ2RGB(XYZ.Gdhp, illuminant = 'E')
HSV.Gdhp <- RGB2HSV(RGB.Gdhp)
hue.Gdhp <- HSV.Gdhp[, 1] * 360
names(hue.Gdhp) <- NULL

spec.Sdhp <-
  subset(
    rev_spectra,
    sample == "Sdhp" &
      variable > 359 & variable < 800,
    select = c(3, 4)
  )
spec.Sdhp <- spec.Sdhp[seq(1, nrow(spec.Sdhp), 5), ]
spec.Sdhp$value <-  (10 ^ -spec.Sdhp$value) * 100
spec.Sdhp <- data.matrix(spec.Sdhp)
XYZ.Sdhp <-
  spectra2XYZ(spec.Sdhp, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.Sdhp <- XYZ2RGB(XYZ.Sdhp, illuminant = 'E')
HSV.Sdhp <- RGB2HSV(RGB.Sdhp)
hue.Sdhp <- HSV.Sdhp[, 1] * 360
names(hue.Sdhp) <- NULL

spec.Hdhp <-
  subset(
    rev_spectra,
    sample == "Hdhp" &
      variable > 359 & variable < 800,
    select = c(3, 4)
  )
spec.Hdhp <- spec.Hdhp[seq(1, nrow(spec.Hdhp), 5), ]
spec.Hdhp$value <-  (10 ^ -spec.Hdhp$value) * 100
spec.Hdhp <- data.matrix(spec.Hdhp)
XYZ.Hdhp <-
  spectra2XYZ(spec.Hdhp, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.Hdhp <- XYZ2RGB(XYZ.Hdhp, illuminant = 'E')
HSV.Hdhp <- RGB2HSV(RGB.Hdhp)
hue.Hdhp <- HSV.Hdhp[, 1] * 360
names(hue.Hdhp) <- NULL