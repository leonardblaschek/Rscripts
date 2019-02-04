library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(reshape2)
library(showtext)
library(colorscience)
library(cowplot)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

liq.DHP <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-12-06_DHPs/17-12-06_DHPs_stained_3min_mixed.csv"
  )
liq.DHP <- melt(liq.DHP, id = c("Well"))
liq.DHP$variable <- gsub("X([0-9]+)", "\\1", liq.DHP$variable)
liq.DHP$variable <- as.numeric(liq.DHP$variable)
liq.DHP[, 4] <- 'stained'
colnames(liq.DHP)[4] <- 'condition'
colnames(liq.DHP)[1] <- 'compound'

liq.DHP.unstained <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-12-06_DHPs/17-12-06_DHPs_unstained.csv"
  )
liq.DHP.unstained <- melt(liq.DHP.unstained, id = c("Well"))
liq.DHP.unstained$variable <-
  gsub("X([0-9]+)", "\\1", liq.DHP.unstained$variable)
liq.DHP.unstained$variable <- as.numeric(liq.DHP.unstained$variable)
liq.DHP.unstained[, 4] <- 'unstained'
colnames(liq.DHP.unstained)[4] <- 'condition'
colnames(liq.DHP.unstained)[1] <- 'compound'

# subtract blank
liq.DHP <- rbind(liq.DHP, liq.DHP.unstained)
liq.DHP[, 1] <-
  revalue(
    liq.DHP[, 1],
    c(
      'A1' = 'G-OH',
      'A2' = 'G-COOH',
      'A3' = 'G-CHO',
      'A4' = 'S-OH',
      'A5' = 'S-COOH',
      'A6' = 'S-CHO',
      'A7' = 'blank'
    )
  )
liq.DHP.blank <-
  subset(liq.DHP, compound == 'blank', select = c(2, 3, 4))
colnames(liq.DHP.blank)[2] <- 'blank'
liq.DHP <-
  merge(liq.DHP,
        liq.DHP.blank,
        by = c("variable", "condition"),
        all = TRUE)
liq.DHP$value <- liq.DHP$value - liq.DHP$blank
liq.DHP[, 5] <- NULL

# calculate HSV colours for G-CHO
data(illuminants)
spec.GCHO <-
  subset(
    liq.DHP,
    compound == 'G-CHO' &
      condition == 'stained' & variable > 359 & variable < 800,
    select = c(1, 4)
  )
spec.GCHO <- spec.GCHO[seq(1, nrow(spec.GCHO), 5),]
spec.GCHO$value <-  (10 ^ -spec.GCHO$value) * 100
spec.GCHO <- data.matrix(spec.GCHO)
XYZ.GCHO <-
  spectra2XYZ(spec.GCHO, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.GCHO <- XYZ2RGB(XYZ.GCHO, illuminant = 'E')
HSV.GCHO <- RGB2HSV(RGB.GCHO)
hue.GCHO <- HSV.GCHO[, 1] * 360
names(hue.GCHO) <- NULL

# calculate HSV colours for S-CHO
spec.SCHO <-
  subset(
    liq.DHP,
    compound == 'S-CHO' &
      condition == 'stained' & variable > 359 & variable < 800,
    select = c(1, 4)
  )
spec.SCHO <- spec.SCHO[seq(1, nrow(spec.SCHO), 5),]
spec.SCHO$value <-  (10 ^ -spec.SCHO$value) * 100
spec.SCHO <- data.matrix(spec.SCHO)
XYZ.SCHO <-
  spectra2XYZ(spec.SCHO, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.SCHO <- XYZ2RGB(XYZ.SCHO, illuminant = 'E')
HSV.SCHO <- RGB2HSV(RGB.SCHO)
hue.SCHO <- HSV.SCHO[, 1] * 360
names(hue.SCHO) <- NULL

# labels for plot
labels <-
  data.frame(280,
             "CHO",
             "stained",
             "G",
             3.2,
             round(hue.GCHO, 0),
             "GCHO")
labels[, 4] <- factor("G", levels = c("G", "S"))
labels[, 7] <- factor("GCHO", levels = c("GCHO", "SCHO"))
labels[2, ] <-
  c(280, "CHO", "stained", "S", 3.2, round(hue.SCHO, 0), "SCHO")
colnames(labels) <-
  c("variable",
    "gamma",
    "condition",
    "ring",
    "value",
    "lbl",
    "fll")
labels$variable <- as.numeric(as.character(labels$variable))
labels$value <- as.numeric(as.character(labels$value))

# divide by 280 nm
liq.DHP.280 <-
  subset(liq.DHP, variable == '280', select = c(2, 3, 4))
colnames(liq.DHP.280)[3] <- 'reference'
liq.DHP <-
  merge(liq.DHP,
        liq.DHP.280,
        # by = c("condition", "compound"),
        all = TRUE)
liq.DHP <- subset(liq.DHP, compound != 'blank')
liq.DHP$value <- liq.DHP$value / liq.DHP$reference
liq.DHP <- liq.DHP %>%
  select(-reference)

liq.DHP$compound <-
  ordered(liq.DHP$compound,
          levels = c("G-CHO", "S-CHO", "G-COOH", "S-COOH", "G-OH", "S-OH"))
liq.DHP$ring[liq.DHP$compound == "G-CHO"] <- "G"
liq.DHP$ring[liq.DHP$compound == "G-COOH"] <- "G"
liq.DHP$ring[liq.DHP$compound == "G-OH"] <- "G"

liq.DHP$ring[liq.DHP$compound == "S-CHO"] <- "S"
liq.DHP$ring[liq.DHP$compound == "S-COOH"] <- "S"
liq.DHP$ring[liq.DHP$compound == "S-OH"] <- "S"

liq.DHP$gamma[liq.DHP$compound == "G-CHO"] <- "CHO"
liq.DHP$gamma[liq.DHP$compound == "S-CHO"] <- "CHO"

liq.DHP$gamma[liq.DHP$compound == "G-COOH"] <- "COOH"
liq.DHP$gamma[liq.DHP$compound == "S-COOH"] <- "COOH"

liq.DHP$gamma[liq.DHP$compound == "G-OH"] <- "OH"
liq.DHP$gamma[liq.DHP$compound == "S-OH"] <- "OH"

g.dhp.liq <- subset(liq.DHP, ring == "G")
colnames(g.dhp.liq)[3] <- "wavelength"
g.dhp.liq$state <- "liquid"

G.DHP <- merge(g.dhp.liq, g.dhp.solid, by = c("wavelength", "condition", "gamma", "state", "value"), all = TRUE)

# DSL <-
#   ggplot(liq.DHP,
#          aes(y = value, x = variable, group = condition)) +
#   # geom_vline(xintercept = 558,
#   #            color = "grey35",
#   #            linetype = 1) +
#   # annotate(
#   #   "text",
#   #   x = 567,
#   #   y = 2.6,
#   #   label = '558',
#   #   size = 3,
#   #   hjust = 0,
#   #   color = "grey35",
# #   family = 'Helvetica'
# # ) +
# geom_rect(
#   data = labels,
#   aes(fill = fll),
#   xmin = -Inf,
#   xmax = Inf,
#   ymin = -Inf,
#   ymax = Inf,
#   alpha = 0.5
# ) +
#   scale_fill_manual(values = c(
#     "GCHO" = hsv(HSV.GCHO[, 1], HSV.GCHO[, 2], 1),
#     "SCHO" = hsv(HSV.SCHO[, 1], HSV.SCHO[, 2] , 1)
#   ), guide = FALSE) +
#   geom_line(aes(linetype = condition)) +
#   #     scale_y_continuous(breaks=c(0,2)) +
#   theme_minimal() +
#   theme(
#     text = element_text(size = 15, family = "Helvetica", colour = "black"),
#     strip.text = element_text(hjust = 0, face = "italic"),
#     # axis.line.y = element_line(size = 0.75, lineend = "square"),
#     axis.ticks = element_line(
#       size = 0.25,
#       lineend = "square",
#       color = "black"
#     ),
#     axis.title = element_text(size = 12),
#     axis.text.y = element_text(size = 10, colour = "black"),
#     axis.text.x = element_text(
#       size = 10,
#       angle = 0,
#       vjust = 1,
#       hjust = 0.5, 
#       colour = "black"
#     ),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_rect(fill = NA, color = "black", size = 0.25),
#     panel.spacing = unit(1.5, "mm"),
#     # plot.margin = unit(c(0, 0, 0, 0), "cm"),
#     legend.position = c(0.35,0.27),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 9, colour = "black"),
#     legend.key.height = unit(4, "mm")
#   ) +
#   xlim(270, 650) +
#   scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(-0.1, 3.5)) +
#   labs(x = "Wavelength [nm]", y = "Normalised absorbance") +
#   geom_text(
#     data = labels,
#     aes(label = paste("calc. hue: ", lbl)),
#     family = "Helvetica",
#     colour = "black",
#     size = 4,
#     hjust = 0
#   ) +
#   facet_grid(gamma ~ ring)
# pdf("spectra_DHPs_normalised.pdf",
#     height = 4,
#     width = 4)
# DSL
# dev.off()

labels <-
  data.frame(280,
             "CHO",
             "stained",
             "liquid",
             1,
             round(hue.GCHO, 0),
             "Gliq")
labels[, 4] <- factor("liquid", levels = c("liquid", "solid"))
labels[, 7] <- factor("Gliq", levels = c("Gliq", "Gsol"))
labels[2, ] <-
  c(280, "CHO", "stained", "solid", 2, round(hue.CHODS, 0), "Gsol")
colnames(labels) <-
  c("wavelength",
    "gamma",
    "condition",
    "state",
    "value",
    "lbl",
    "fll")
labels$wavelength <- as.numeric(as.character(labels$wavelength))
labels$value <- as.numeric(as.character(labels$value))

DSLS <-
  ggplot(data = G.DHP,
         aes(y = value, x = wavelength, group = condition)) +
# geom_vline(xintercept = 558,
           # color = "grey35",
           # linetype = 1) +
# annotate(
#   "text",
#   x = 567,
#   y = 2.6,
#   label = '558',
#   size = 3,
#   hjust = 0,
#   color = "grey35",
#   family = 'Helvetica'
# ) +
geom_rect(
  data = labels,
  aes(fill = fll),
  xmin = -Inf,
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.5
) +
scale_fill_manual(values = c(
  "Gliq" = hsv(HSV.GCHO[, 1], 0.75, 1),
  "Gsol" = hsv(HSV.CHODS[, 1], 0.75, 1)
), guide = FALSE) +
  geom_line(aes(linetype = condition), size = 0.25) +
  #     scale_y_continuous(breaks=c(0,2)) +
  theme_minimal() +
  theme(
    text = element_text(size = 15, family = "Helvetica", colour = "black"),
    strip.text = element_text(hjust = 0, face = "italic"),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.x = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title = element_text(size = 12),
    # axis.text.y = element_blank(),
    # axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(
      size = 10,
      angle = 0,
      vjust = 1,
      hjust = 0.5, 
      colour = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing = unit(1.5, "mm"),
    # plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, colour = "black"),
    legend.key.height = unit(4, "mm")
  ) +
  xlim(270, 650) +
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-0.1, 2.6)) +
  labs(x = "Wavelength [nm]", y = "Absorbance") +
  # geom_text(
  #   data = labels,
  #   aes(label = paste("calc. hue: ", lbl)),
  #   family = "Helvetica",
  #   colour = "black",
  #   size = 4,
  #   hjust = 0
  # ) +
  facet_grid(state ~ gamma,
             # scales = "free_y"
  )
pdf("spectra_DHPs.pdf",
    height = 3.75,
    width = 5)
DSLS
dev.off()
