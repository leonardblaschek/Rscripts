library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(reshape2)
library(showtext)
library(colorscience)
library(zoo)

#### import Helvetica Neue ####
font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

#### import data for G-DHPs ####

GOHd_UN00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DOHU00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHd_UN00$absorption <-
  GOHd_UN00$absorption + (0 - subset(GOHd_UN00, wavelength ==
                                       700)$absorption)
GOHd_UN01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DOHU01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHd_UN01$absorption <-
  GOHd_UN01$absorption + (0 - subset(GOHd_UN01, wavelength ==
                                       700)$absorption)
GOHd_UN02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DOHU02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHd_UN02$absorption <-
  GOHd_UN02$absorption + (0 - subset(GOHd_UN02, wavelength ==
                                       700)$absorption)

GOHd_ST00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DOHS00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHd_ST00$absorption <-
  GOHd_ST00$absorption + (0 - subset(GOHd_ST00, wavelength ==
                                       700)$absorption)
GOHd_ST01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DOHS01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHd_ST01$absorption <-
  GOHd_ST01$absorption + (0 - subset(GOHd_ST01, wavelength ==
                                       700)$absorption)
GOHd_ST02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DOHS02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHd_ST02$absorption <-
  GOHd_ST02$absorption + (0 - subset(GOHd_ST02, wavelength ==
                                       700)$absorption)

GCHOd_UN00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCHOU00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOd_UN00$absorption <-
  GCHOd_UN00$absorption + (0 - subset(GCHOd_UN00, wavelength ==
                                        700)$absorption)
GCHOd_UN01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCHOU01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOd_UN01$absorption <-
  GCHOd_UN01$absorption + (0 - subset(GCHOd_UN01, wavelength ==
                                        700)$absorption)
GCHOd_UN02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCHOU02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOd_UN02$absorption <-
  GCHOd_UN02$absorption + (0 - subset(GCHOd_UN02, wavelength ==
                                        700)$absorption)

GCHOd_ST00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCHOS00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOd_ST00$absorption <-
  GCHOd_ST00$absorption + (0 - subset(GCHOd_ST00, wavelength ==
                                        700)$absorption)
GCHOd_ST01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCHOS01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOd_ST01$absorption <-
  GCHOd_ST01$absorption + (0 - subset(GCHOd_ST01, wavelength ==
                                        700)$absorption)
GCHOd_ST02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCHOS02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOd_ST02$absorption <-
  GCHOd_ST02$absorption + (0 - subset(GCHOd_ST02, wavelength ==
                                        700)$absorption)

GCOOHd_UN00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHU00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHd_UN00$absorption <-
  GCOOHd_UN00$absorption + (0 - subset(GCOOHd_UN00, wavelength ==
                                         700)$absorption)
GCOOHd_UN01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHU01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHd_UN01$absorption <-
  GCOOHd_UN01$absorption + (0 - subset(GCOOHd_UN01, wavelength ==
                                         700)$absorption)
GCOOHd_UN02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHU02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHd_UN02$absorption <-
  GCOOHd_UN02$absorption + (0 - subset(GCOOHd_UN02, wavelength ==
                                         700)$absorption)

GCOOHd_ST00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHS00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHd_ST00$absorption <-
  GCOOHd_ST00$absorption + (0 - subset(GCOOHd_ST00, wavelength ==
                                         700)$absorption)
GCOOHd_ST01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHS01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHd_ST01$absorption <-
  GCOOHd_ST01$absorption + (0 - subset(GCOOHd_ST01, wavelength ==
                                         700)$absorption)
GCOOHd_ST02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHS02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHd_ST02$absorption <-
  GCOOHd_ST02$absorption + (0 - subset(GCOOHd_ST02, wavelength ==
                                         700)$absorption)

#### unify data ####
GOHd_UN <- rbind(GOHd_UN00, GOHd_UN01, GOHd_UN02)
GOHd_UN$residue <- "G-OH"
GOHd_UN$structure <- "DHP"
GOHd_UN$condition <- "unstained"
GOHd_ST <- rbind(GOHd_ST00, GOHd_ST01, GOHd_ST02)
GOHd_ST$residue <- "G-OH"
GOHd_ST$structure <- "DHP"
GOHd_ST$condition <- "stained"
GOHd_ST$absorption[GOHd_ST$absorption > 4.9] <- NA
GCHOd_UN <- rbind(GCHOd_UN00, GCHOd_UN01, GCHOd_UN02)
GCHOd_UN$residue <- "G-CHO"
GCHOd_UN$structure <- "DHP"
GCHOd_UN$condition <- "unstained"
GCHOd_UN$absorption[GCHOd_UN$absorption > 4.9] <- NA
GCHOd_ST <- rbind(GCHOd_ST00, GCHOd_ST01, GCHOd_ST02)
GCHOd_ST$residue <- "G-CHO"
GCHOd_ST$structure <- "DHP"
GCHOd_ST$condition <- "stained"
GCHOd_ST$absorption[GCHOd_ST$absorption > 4.9] <- NA
GCOOHd_UN <- rbind(GCOOHd_UN00, GCOOHd_UN01, GCOOHd_UN02)
GCOOHd_UN$residue <- "G-COOH"
GCOOHd_UN$structure <- "DHP"
GCOOHd_UN$condition <- "unstained"
GCOOHd_ST <- rbind(GCOOHd_ST00, GCOOHd_ST01, GCOOHd_ST02)
GCOOHd_ST$residue <- "G-COOH"
GCOOHd_ST$structure <- "DHP"
GCOOHd_ST$condition <- "stained"
GCOOHd_ST$absorption[GCOOHd_ST$absorption > 4.9] <- NA

spectra.Gd <-
  ddply(
    GOHd_UN,
    c("wavelength", "residue", "condition", "structure"),
    summarise,
    mean = mean(absorption)
  )
spectra.Gd <-
  rbind(spectra.Gd, (ddply(
    GOHd_ST,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption, na.rm = TRUE)
  )))
spectra.Gd <-
  rbind(spectra.Gd, (ddply(
    GCHOd_UN,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption, na.rm = TRUE)
  )))
spectra.Gd <-
  rbind(spectra.Gd, (ddply(
    GCHOd_ST,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption, na.rm = TRUE)
  )))
spectra.Gd <-
  rbind(spectra.Gd, (ddply(
    GCOOHd_UN,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption, na.rm = TRUE)
  )))
spectra.Gd <-
  rbind(spectra.Gd, (ddply(
    GCOOHd_ST,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption, na.rm = TRUE)
  )))

spectra <- spectra.Gd
spectra$time <- '0 h'
spectra$smooth <-
  rollapply(
    spectra$mean,
    3,
    mean,
    fill = "extend",
    align = "center",
    partial = TRUE
  )

write.csv(spectra, "DHP_spectra.csv")

#### Colour Calculations ####
data(illuminants)

spec.CHODS <-
  subset(
    spectra,
    structure == 'DHP' &
      residue == 'G-CHO' &
      condition == 'stained' & wavelength > 359,
    select = c(1, 5)
  )
spec.CHODS <- spec.CHODS[seq(1, nrow(spec.CHODS), 5), ]
spec.CHODS$mean <-  (10 ^ -spec.CHODS$mean) * 100
spec.CHODS <- data.matrix(spec.CHODS)
XYZ.CHODS <-
  spectra2XYZ(spec.CHODS, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.CHODS <- XYZ2RGB(XYZ.CHODS, illuminant = 'E')
HSV.CHODS <- RGB2HSV(RGB.CHODS)
hue.CHODS <- HSV.CHODS[, 1] * 360
names(hue.CHODS) <- NULL

#### prepare solid spectra df for merging with liquid spectra df ####
g.dhp.solid <- subset(spectra, structure == "DHP")
g.dhp.solid$gamma[g.dhp.solid$residue == "G-CHO"] <- "CHO"
g.dhp.solid$gamma[g.dhp.solid$residue == "G-COOH"] <- "COOH"
g.dhp.solid$gamma[g.dhp.solid$residue == "G-OH"] <- "OH"
g.dhp.solid$state <- "solid"
colnames(g.dhp.solid)[7] <- "value"

#### read liquid spectra data ####
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

#### subtract blank ####
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

#### colour calculations ####
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

liq.DHP <- subset(liq.DHP, compound != 'blank')

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
colnames(g.dhp.liq)[1] <- "wavelength"
g.dhp.liq$state <- "liquid"

#### merge liquid and solid spectra dfs ####
G.DHP <-
  merge(
    g.dhp.liq,
    g.dhp.solid,
    by = c("wavelength", "condition", "gamma", "state", "value"),
    all = TRUE
  )

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

#### plot spectra ####
DSLS <-
  ggplot(data = G.DHP,
         aes(y = value, x = wavelength, group = condition)) +
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
  ),
  guide = FALSE) +
  geom_line(aes(linetype = condition), size = 0.25) +
  theme_minimal() +
  theme(
    text = element_text(
      size = 15,
      family = "Helvetica",
      colour = "black"
    ),
    strip.text = element_text(hjust = 0, face = "italic"),
    axis.ticks.x = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title = element_text(size = 12),
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
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, colour = "black"),
    legend.key.height = unit(4, "mm")
  ) +
  xlim(270, 650) +
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-0.1, 2.5)) +
  labs(x = "Wavelength [nm]", y = "Absorbance") +
  facet_grid(state ~ gamma)
pdf("spectra_DHPs.pdf",
    height = 3.75,
    width = 5)
DSLS
dev.off()