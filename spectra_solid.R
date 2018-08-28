library(ggplot2)
library(ggthemes)
library(plyr)
library(reshape2)
library(graphics)
library(zoo)
library(colorscience)

# ---- G-monomers ----

GOHm_UN00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/OH_UN00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHm_UN00$absorption <-
  GOHm_UN00$absorption + (0 - subset(GOHm_UN00, wavelength ==
                                       700)$absorption)
GOHm_UN01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/OH_UN01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHm_UN01$absorption <-
  GOHm_UN01$absorption + (0 - subset(GOHm_UN01, wavelength ==
                                       700)$absorption)
GOHm_UN02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/OH_UN02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHm_UN02$absorption <-
  GOHm_UN02$absorption + (0 - subset(GOHm_UN02, wavelength ==
                                       700)$absorption)

GOHm_ST00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/OH_ST00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHm_ST00$absorption <-
  GOHm_ST00$absorption + (0 - subset(GOHm_ST00, wavelength ==
                                       700)$absorption)
GOHm_ST01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/OH_ST01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHm_ST01$absorption <-
  GOHm_ST01$absorption + (0 - subset(GOHm_ST01, wavelength ==
                                       700)$absorption)
GOHm_ST02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/OH_ST02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GOHm_ST02$absorption <-
  GOHm_ST02$absorption + (0 - subset(GOHm_ST02, wavelength ==
                                       700)$absorption)

GCHOm_UN00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/G'600U00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOm_UN00$absorption <-
  GCHOm_UN00$absorption + (0 - subset(GCHOm_UN00, wavelength ==
                                        700)$absorption)
GCHOm_UN01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/G'600U01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOm_UN01$absorption <-
  GCHOm_UN01$absorption + (0 - subset(GCHOm_UN01, wavelength ==
                                        700)$absorption)
GCHOm_UN02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/G'600U02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOm_UN02$absorption <-
  GCHOm_UN02$absorption + (0 - subset(GCHOm_UN02, wavelength ==
                                        700)$absorption)

GCHOm_ST00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/G'600S00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOm_ST00$absorption <-
  GCHOm_ST00$absorption + (0 - subset(GCHOm_ST00, wavelength ==
                                        700)$absorption)
GCHOm_ST01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/G'600S01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOm_ST01$absorption <-
  GCHOm_ST01$absorption + (0 - subset(GCHOm_ST01, wavelength ==
                                        700)$absorption)
GCHOm_ST02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/G'600S02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCHOm_ST02$absorption <-
  GCHOm_ST02$absorption + (0 - subset(GCHOm_ST02, wavelength ==
                                        700)$absorption)

GCOOHm_UN00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/COOH_U00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHm_UN00$absorption <-
  GCOOHm_UN00$absorption + (0 - subset(GCOOHm_UN00, wavelength ==
                                         700)$absorption)
GCOOHm_UN01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/COOH_U01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHm_UN01$absorption <-
  GCOOHm_UN01$absorption + (0 - subset(GCOOHm_UN01, wavelength ==
                                         700)$absorption)
GCOOHm_UN02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/COOH_U02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHm_UN02$absorption <-
  GCOOHm_UN02$absorption + (0 - subset(GCOOHm_UN02, wavelength ==
                                         700)$absorption)

GCOOHm_ST00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/COOH_S00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHm_ST00$absorption <-
  GCOOHm_ST00$absorption + (0 - subset(GCOOHm_ST00, wavelength ==
                                         700)$absorption)
GCOOHm_ST01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/COOH_S01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHm_ST01$absorption <-
  GCOOHm_ST01$absorption + (0 - subset(GCOOHm_ST01, wavelength ==
                                         700)$absorption)
GCOOHm_ST02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/600muM/COOH_S02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
GCOOHm_ST02$absorption <-
  GCOOHm_ST02$absorption + (0 - subset(GCOOHm_ST02, wavelength ==
                                         700)$absorption)


# GOHm_ST00$absorption <- GOHm_ST00$absorption + (1 - subset(GOHm_UN00, wavelength ==
#     280)$absorption)
# GOHm_ST01$absorption <- GOHm_ST01$absorption + (1 - subset(GOHm_UN01, wavelength ==
#     280)$absorption)
# GOHm_ST02$absorption <- GOHm_ST02$absorption + (1 - subset(GOHm_UN02, wavelength ==
#     280)$absorption)
#
# GOHm_UN00$absorption <- GOHm_UN00$absorption + (1 - subset(GOHm_UN00, wavelength ==
#     280)$absorption)
# GOHm_UN01$absorption <- GOHm_UN01$absorption + (1 - subset(GOHm_UN01, wavelength ==
#     280)$absorption)
# GOHm_UN02$absorption <- GOHm_UN02$absorption + (1 - subset(GOHm_UN02, wavelength ==
#     280)$absorption)
#
# GCHOm_ST00$absorption <- GCHOm_ST00$absorption + (1 - subset(GCHOm_UN00, wavelength ==
#     280)$absorption)
# GCHOm_ST01$absorption <- GCHOm_ST01$absorption + (1 - subset(GCHOm_UN01, wavelength ==
#     280)$absorption)
# GCHOm_ST02$absorption <- GCHOm_ST02$absorption + (1 - subset(GCHOm_UN02, wavelength ==
#     280)$absorption)
#
# GCHOm_UN00$absorption <- GCHOm_UN00$absorption + (1 - subset(GCHOm_UN00, wavelength ==
#     280)$absorption)
# GCHOm_UN01$absorption <- GCHOm_UN01$absorption + (1 - subset(GCHOm_UN01, wavelength ==
#     280)$absorption)
# GCHOm_UN02$absorption <- GCHOm_UN02$absorption + (1 - subset(GCHOm_UN02, wavelength ==
#     280)$absorption)
#
# GCOOHm_ST00$absorption <- GCOOHm_ST00$absorption + (1 - subset(GCOOHm_UN00, wavelength ==
#     280)$absorption)
# GCOOHm_ST01$absorption <- GCOOHm_ST01$absorption + (1 - subset(GCOOHm_UN01, wavelength ==
#     280)$absorption)
# GCOOHm_ST02$absorption <- GCOOHm_ST02$absorption + (1 - subset(GCOOHm_UN02, wavelength ==
#     280)$absorption)
#
# GCOOHm_UN00$absorption <- GCOOHm_UN00$absorption + (1 - subset(GCOOHm_UN00, wavelength ==
#     280)$absorption)
# GCOOHm_UN01$absorption <- GCOOHm_UN01$absorption + (1 - subset(GCOOHm_UN01, wavelength ==
#     280)$absorption)
# GCOOHm_UN02$absorption <- GCOOHm_UN02$absorption + (1 - subset(GCOOHm_UN02, wavelength ==
#     280)$absorption)

# ---- G-DHPs ----

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


# GOHd_ST00$absorption <- GOHd_ST00$absorption + (1 - subset(GOHd_UN00, wavelength ==
#     280)$absorption)
# GOHd_ST01$absorption <- GOHd_ST01$absorption + (1 - subset(GOHd_UN01, wavelength ==
#     280)$absorption)
# GOHd_ST02$absorption <- GOHd_ST02$absorption + (1 - subset(GOHd_UN02, wavelength ==
#     280)$absorption)
#
# GOHd_UN00$absorption <- GOHd_UN00$absorption + (1 - subset(GOHd_UN00, wavelength ==
#     280)$absorption)
# GOHd_UN01$absorption <- GOHd_UN01$absorption + (1 - subset(GOHd_UN01, wavelength ==
#     280)$absorption)
# GOHd_UN02$absorption <- GOHd_UN02$absorption + (1 - subset(GOHd_UN02, wavelength ==
#     280)$absorption)
#
# GCHOd_ST00$absorption <- GCHOd_ST00$absorption + (1 - subset(GCHOd_UN00, wavelength ==
#     280)$absorption)
# GCHOd_ST01$absorption <- GCHOd_ST01$absorption + (1 - subset(GCHOd_UN01, wavelength ==
#     280)$absorption)
# GCHOd_ST02$absorption <- GCHOd_ST02$absorption + (1 - subset(GCHOd_UN02, wavelength ==
#     280)$absorption)
#
# GCHOd_UN00$absorption <- GCHOd_UN00$absorption + (1 - subset(GCHOd_UN00, wavelength ==
#     280)$absorption)
# GCHOd_UN01$absorption <- GCHOd_UN01$absorption + (1 - subset(GCHOd_UN01, wavelength ==
#     280)$absorption)
# GCHOd_UN02$absorption <- GCHOd_UN02$absorption + (1 - subset(GCHOd_UN02, wavelength ==
#     280)$absorption)
#
# GCOOHd_ST00$absorption <- GCOOHd_ST00$absorption + (1 - subset(GCOOHd_UN00, wavelength ==
#     280)$absorption)
# GCOOHd_ST01$absorption <- GCOOHd_ST01$absorption + (1 - subset(GCOOHd_UN01, wavelength ==
#     280)$absorption)
# GCOOHd_ST02$absorption <- GCOOHd_ST02$absorption + (1 - subset(GCOOHd_UN02, wavelength ==
#     280)$absorption)
#
# GCOOHd_UN00$absorption <- GCOOHd_UN00$absorption + (1 - subset(GCOOHd_UN00, wavelength ==
#     280)$absorption)
# GCOOHd_UN01$absorption <- GCOOHd_UN01$absorption + (1 - subset(GCOOHd_UN01, wavelength ==
#     280)$absorption)
# GCOOHd_UN02$absorption <- GCOOHd_UN02$absorption + (1 - subset(GCOOHd_UN02, wavelength ==
#     280)$absorption)

# normalise <- function(x, y) {x$absorption <- x$absorption/y$absorption[y$wavelength == 280]}
# 
# GOHm_ST00$absorption <- normalise(GOHm_ST00, GOHm_UN00)
# GOHm_ST01$absorption <- normalise(GOHm_ST01, GOHm_UN01)
# GOHm_ST02$absorption <- normalise(GOHm_ST02, GOHm_UN02)
# 
# GCOOHm_ST00$absorption <- normalise(GCOOHm_ST00, GCOOHm_UN00)
# GCOOHm_ST01$absorption <- normalise(GCOOHm_ST01, GCOOHm_UN01)
# GCOOHm_ST02$absorption <- normalise(GCOOHm_ST02, GCOOHm_UN02)
# 
# GCHOm_ST00$absorption <- normalise(GCHOm_ST00, GCHOm_UN00)
# GCHOm_ST01$absorption <- normalise(GCHOm_ST01, GCHOm_UN01)
# GCHOm_ST02$absorption <- normalise(GCHOm_ST02, GCHOm_UN02)
# 
# GOHd_ST00$absorption <- normalise(GOHd_ST00, GOHd_UN00)
# GOHd_ST01$absorption <- normalise(GOHd_ST01, GOHd_UN01)
# GOHd_ST02$absorption <- normalise(GOHd_ST02, GOHd_UN02)
# 
# GCOOHd_ST00$absorption <- normalise(GCOOHd_ST00, GCOOHd_UN00)
# GCOOHd_ST01$absorption <- normalise(GCOOHd_ST01, GCOOHd_UN01)
# GCOOHd_ST02$absorption <- normalise(GCOOHd_ST02, GCOOHd_UN02)
# 
# GCHOd_ST00$absorption <- normalise(GCHOd_ST00, GCHOd_UN00)
# GCHOd_ST01$absorption <- normalise(GCHOd_ST01, GCHOd_UN01)
# GCHOd_ST02$absorption <- normalise(GCHOd_ST02, GCHOd_UN02)

GOHm_UN <- rbind(GOHm_UN00, GOHm_UN01, GOHm_UN02)
GOHm_UN$residue <- "G-OH"
GOHm_UN$structure <- "Monomer"
GOHm_UN$condition <- "unstained"
GOHm_ST <- rbind(GOHm_ST00, GOHm_ST01, GOHm_ST02)
GOHm_ST$residue <- "G-OH"
GOHm_ST$structure <- "Monomer"
GOHm_ST$condition <- "stained"
GCHOm_UN <- rbind(GCHOm_UN00, GCHOm_UN01, GCHOm_UN02)
GCHOm_UN$residue <- "G-CHO"
GCHOm_UN$structure <- "Monomer"
GCHOm_UN$condition <- "unstained"
GCHOm_ST <- rbind(GCHOm_ST00, GCHOm_ST01, GCHOm_ST02)
GCHOm_ST$residue <- "G-CHO"
GCHOm_ST$structure <- "Monomer"
GCHOm_ST$condition <- "stained"
GCOOHm_UN <- rbind(GCOOHm_UN00, GCOOHm_UN01, GCOOHm_UN02)
GCOOHm_UN$residue <- "G-COOH"
GCOOHm_UN$structure <- "Monomer"
GCOOHm_UN$condition <- "unstained"
GCOOHm_ST <- rbind(GCOOHm_ST00, GCOOHm_ST01, GCOOHm_ST02)
GCOOHm_ST$residue <- "G-COOH"
GCOOHm_ST$structure <- "Monomer"
GCOOHm_ST$condition <- "stained"
spectra.Gm <-
  ddply(
    GOHm_UN,
    c("wavelength", "residue", "condition", "structure"),
    summarise,
    mean = mean(absorption)
  )
spectra.Gm <-
  rbind(spectra.Gm, (ddply(
    GOHm_ST,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption)
  )))
spectra.Gm <-
  rbind(spectra.Gm, (ddply(
    GCHOm_UN,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption)
  )))
spectra.Gm <-
  rbind(spectra.Gm, (ddply(
    GCHOm_ST,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption)
  )))
spectra.Gm <-
  rbind(spectra.Gm, (ddply(
    GCOOHm_UN,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption)
  )))
spectra.Gm <-
  rbind(spectra.Gm, (ddply(
    GCOOHm_ST,
    c("wavelength", "residue", "condition",
      "structure"),
    summarise,
    mean = mean(absorption)
  )))

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


# ---- Phloroglucinol ----

PHLOG_ST00 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHS00.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
PHLOG_ST00$absorption <-
  PHLOG_ST00$absorption + (0 - subset(PHLOG_ST00, wavelength ==
                                        700)$absorption)
PHLOG_ST01 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHS01.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
PHLOG_ST01$absorption <-
  PHLOG_ST01$absorption + (0 - subset(PHLOG_ST01, wavelength ==
                                        700)$absorption)
PHLOG_ST02 <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_DHP_spectra/DHP/DCOOHS02.ASC",
    col.names = c("wavelength", "absorption"),
    header = FALSE
  )
PHLOG_ST02$absorption <-
  PHLOG_ST02$absorption + (0 - subset(PHLOG_ST02, wavelength ==
                                        700)$absorption)

PHLOG_ST <- rbind(PHLOG_ST00, PHLOG_ST01, PHLOG_ST02)
PHLOG_ST$residue <- "Phloroglucinol-HCl"
PHLOG_ST$structure <- "Monomer"
PHLOG_ST$condition <- "stained"

spectra <- rbind(spectra.Gm, spectra.Gd)
spectra$time <- '0 h'
# spectra <- rbind(spectra, spectra.fade)
# spectra <- rbind(spectra,(ddply(PHLOG_ST, c("wavelength", "residue", "condition", "structure"), summarise, mean = mean(absorption, na.rm=TRUE))))

# GCHOm.342 <- spectra$mean[spectra$residue=="G-CHO" & spectra$condition=="unstained" &
#      spectra$wavelength==342 & spectra$structure=="Monomer"]
# GCHOd.342 <- spectra$mean[spectra$residue=="G-CHO" & spectra$condition=="unstained" &
#      spectra$wavelength==342 & spectra$structure=="DHP"]
# spectra$mean <- ifelse((spectra$residue == "G-CHO" & spectra$condition=="stained" & spectra$structure=="Monomer"), spectra$mean / GCHOm.342, spectra$mean)
# spectra$mean <- ifelse((spectra$residue == "G-CHO" & spectra$condition=="stained" & spectra$structure=="DHP"), spectra$mean / GCHOd.342, spectra$mean)
#
# GOHm.342 <- spectra$mean[spectra$residue=="G-OH" & spectra$condition=="unstained" &
#      spectra$wavelength==342 & spectra$structure=="Monomer"]
# GOHd.342 <- spectra$mean[spectra$residue=="G-OH" & spectra$condition=="unstained" &
#      spectra$wavelength==342 & spectra$structure=="DHP"]
# spectra$mean <- ifelse((spectra$residue == "G-OH" & spectra$condition=="stained" & spectra$structure=="Monomer"), spectra$mean / GOHm.342, spectra$mean)
# spectra$mean <- ifelse((spectra$residue == "G-OH" & spectra$condition=="stained" & spectra$structure=="DHP"), spectra$mean / GOHd.342, spectra$mean)
#
# GCOOHm.342 <- spectra$mean[spectra$residue=="G-COOH" & spectra$condition=="unstained" &
#      spectra$wavelength==342 & spectra$structure=="Monomer"]
# GCOOHd.342 <- spectra$mean[spectra$residue=="G-COOH" & spectra$condition=="unstained" &
#      spectra$wavelength==342 & spectra$structure=="DHP"]
# spectra$mean <- ifelse((spectra$residue == "G-COOH" & spectra$condition=="stained" & spectra$structure=="Monomer"), spectra$mean / GCOOHm.342, spectra$mean)
# spectra$mean <- ifelse((spectra$residue == "G-COOH" & spectra$condition=="stained" & spectra$structure=="DHP"), spectra$mean / GCOOHd.342, spectra$mean)

spectra$smooth <-
  rollapply(
    spectra$mean,
    3,
    mean,
    fill = "extend",
    align = "center",
    partial = TRUE
  )

# ---- Colour Calculations ----

data(illuminants)
spec.CHOMS <-
  subset(
    spectra,
    structure == 'Monomer' &
      residue == 'G-CHO' &
      condition == 'stained' & wavelength > 359,
    select = c(1, 5)
  )
spec.CHOMS <- spec.CHOMS[seq(1, nrow(spec.CHOMS), 5), ]
spec.CHOMS$mean <-  (10 ^ -spec.CHOMS$mean) * 100
spec.CHOMS <- data.matrix(spec.CHOMS)
XYZ.CHOMS <-
  spectra2XYZ(spec.CHOMS, illuminantIn = illuminants[, c('wlnm', 'E')])
RGB.CHOMS <- XYZ2RGB(XYZ.CHOMS, illuminant = 'E')
HSV.CHOMS <- RGB2HSV(RGB.CHOMS)
hue.CHOMS <- HSV.CHOMS[, 1] * 360
names(hue.CHOMS) <- NULL

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

labels <-
  data.frame(280,
             "G-CHO",
             "stained",
             "Monomer",
             2.3,
             "0 h",
             2.3,
             round(hue.CHOMS, 0),
             "CHOMS")
labels[, 4] <- factor("Monomer", levels = c("Monomer", "DHP"))
labels[, 9] <- factor("CHOMS", levels = c("CHOMS", "CHODS"))
labels[2,] <-
  c(280, "G-CHO", "stained", "DHP", 2.3, "0 h", 2.3, round(hue.CHODS, 0), "CHODS")
colnames(labels) <-
  c("wavelength",
    "residue",
    "condition",
    "structure",
    "mean",
    "time",
    "smooth",
    "lbl",
    "fll")
labels$wavelength <- as.numeric(as.character(labels$wavelength))
labels$mean <- as.numeric(as.character(labels$mean))
labels$smooth <- as.numeric(as.character(labels$smooth))

# ---- PLOT ----
# pdf("test.pdf")
# p <- ggplot(GOHd_ST, aes(y=absorption, x=wavelength)) + geom_point()
# p
#dev.off()

g.dhp.solid <- subset(spectra, structure == "DHP")
g.dhp.solid$gamma[g.dhp.solid$residue == "G-CHO"] <- "CHO"
g.dhp.solid$gamma[g.dhp.solid$residue == "G-COOH"] <- "COOH"
g.dhp.solid$gamma[g.dhp.solid$residue == "G-OH"] <- "OH"
g.dhp.solid$state <- "solid"
colnames(g.dhp.solid)[7] <- "value"

DSS <-
  ggplot(
    subset(spectra, time == '0 h'),
    aes(y = smooth, x = wavelength, linetype = condition)
  ) +
  geom_rect(
    data = labels,
    aes(fill = fll),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("CHOMS" = hsv(HSV.CHOMS[,1], 1, 1), "CHODS" = hsv(HSV.CHODS[,1], 1 , 1))) +
  # geom_vline(xintercept = 508,
  #            color = "grey35",
  #            linetype = 1) +
  # geom_vline(xintercept = 560,
  #            color = "grey35",
  #            linetype = 1) +
  #     geom_vline(xintercept = 342, color = "brown", linetype = 1) +
  #     geom_vline(xintercept = 270, color = "pink", linetype = 3) +
  #     geom_hline(yintercept = 0, color="grey35", lwd=0.2) +
  geom_line() +
  theme_minimal() +
  theme(
    text = element_text(size = 15, family = "Helvetica"),
    strip.text.x = element_text(hjust = 0, face = "italic"),
    strip.text.y = element_text(vjust = 0.5, face = "italic"),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10, colour = "black"),
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
    legend.position = "none"
  ) +
  xlim(270, 675) +
  ylim(-.5, 3) +
  scale_color_brewer(palette = "Dark2") +
  geom_text(
    data = labels,
    aes(label = paste("calc. hue: ", lbl)),
    family = "Helvetica",
    colour = "black",
    size = 4,
    hjust = 0
  ) +
      geom_vline(xintercept = 559, size=1) +
  annotate("text", x = 565, y = 2.3, label = '560', size = 3, hjust = 0, color = "grey35", family = 'Helvetica') +
      geom_vline(xintercept = 510, size=1) +
  annotate("text", x = 500, y = 2.3, label = '508', size = 3, hjust = 1, color = "grey35", family = 'Helvetica') +
  labs(x = "Wavelength [nm]", y = "Absorbance") +
  facet_grid(residue ~ structure)

pdf("spectra.pdf", height = 4, width = 4)
DSS
dev.off()
