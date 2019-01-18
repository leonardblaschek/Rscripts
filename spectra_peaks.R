library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(showtext)
library(colorscience)
library(cowplot)
data(illuminants)

#### import Helvetica Neue ####
font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

#### load and tidy data ####
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
peaks.lbl <- peaks[1:12, ]
peaks.lbl$variable <- 280
peaks.lbl$value <- 0.3
peaks <-
  subset(
    peaks,
    peaks$peak == "Peak 1" |
      peaks$peak == "Peak 5" |
      peaks$peak == "Peak 6" |
      peaks$peak == "Peak 9" | peaks$peak == "Peak 12"
  )
peaks.lbl <-
  subset(
    peaks.lbl,
    peaks.lbl$peak == "Peak 1" |
      peaks.lbl$peak == "Peak 5" |
      peaks.lbl$peak == "Peak 6" |
      peaks.lbl$peak == "Peak 9" | peaks.lbl$peak == "Peak 12"
  )

#### calculate HSV colours for Peak 6/9 ####
spec.peak9 <-
  subset(
    peaks,
    peak == "Peak 9" &
      repeat. == 1 &
      condition == "acidified" &
      variable > 359 & variable < 800,
    select = c(4, 5)
  )
spec.peak9 <- spec.peak9[seq(1, nrow(spec.peak9), 5),]
spec.peak9$value <-  (10 ^ -spec.peak9$value) * 100
spec.peak9 <- data.matrix(spec.peak9)
XYZ.peak9 <-
  spectra2XYZ(spec.peak9, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.peak9 <- XYZ2RGB(XYZ.peak9, illuminant = 'E')
HSV.peak9 <- RGB2HSV(RGB.peak9)
hue.peak9 <- HSV.peak9[, 1] * 360
names(hue.peak9) <- NULL

peak9.col <-
  data.frame("Peak 9", "acidified", 1, "p9", 280, 0.225, round(hue.peak9, 0))
colnames(peak9.col) <-
  c("peak",
    "condition",
    "repeat.",
    "fll",
    "variable",
    "value",
    "lbl")

#### plot peak spectra ####
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
  scale_fill_manual(values = c("p9" = hsv(HSV.peak9[, 1], 10 * HSV.peak9[, 2], 1)),
                    guide = FALSE) +
  geom_line() +
  theme_minimal() +
  theme(
    text = element_text(size = 15, family = "Helvetica"),
    strip.text = element_blank(),
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
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.height = unit(4, "mm")
  ) +
  xlim(270, 650) +
  scale_y_continuous(
    breaks = c(0, 0.2),
    limits = c(0, 0.35),
    labels = c("0", "0.2")
  ) +
  labs(x = "Wavelength [nm]", y = "Absorbance") +
  facet_wrap( ~ peak, ncol = 5) +
  geom_text(
    data = peaks.lbl,
    aes(label = peak),
    family = "Helvetica",
    colour = "black",
    fontface = "italic",
    size = 3,
    hjust = 0
  ) +
  geom_text(
    data = peak9.col,
    aes(label = paste("calc. hue: ", lbl)),
    family = "Helvetica",
    colour = "black",
    size = 3,
    hjust = 0
  )

pdf("spectra_peaks.pdf", width = 6, height = 1.5)
SP
dev.off()