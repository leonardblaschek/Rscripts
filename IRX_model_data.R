library(baseline)
library(agricolae)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(showtext)
library(purrr)
library(tibble)
library(zoo)
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

#### generating plot theme ####
theme_leo <- function(base_size = 12,
                      base_family = "Helvetica"){
  theme_minimal(base_size = base_size,
                base_family = base_family) %+replace%
    theme(
      strip.text = element_text(hjust = 0, face = "italic"),
      axis.ticks = element_line(
        size = 0.25,
        lineend = "square",
        color = "black"
      ),
      axis.title = element_blank(),
      axis.text.y = element_text(colour = "black"),
      axis.text.x = element_text(
        colour = "black",
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.25),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "bottom",
      legend.text = element_text(size = rel(0.8)),
      legend.key.height = unit(4, "mm"),
      legend.key.width = unit(30, "mm"),
      # plot.margin = unit(c(0, 0, 0, 0), "cm"),   
      
      complete = TRUE
    )
}

#### functions for scaling and statistics ####
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

tukey <- function(x) {
  aov1 <- aov(data = x, value ~ genotype)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups$groups[["genotype"]] <- rownames(groups$groups)
  groups$groups[["value"]] <- 0
  return(groups[["groups"]])
}

#### import and tidy Wiesner test data ####
wiesner.data <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
           skip = 2)

wiesner.data$genotype <- recode(wiesner.data$genotype, cad4x5 = "cad4xcad5")

# set cell types according to measurement order
wiesner.data[1:50 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "IF"
wiesner.data[51:100 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "MX"
wiesner.data[101:150 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "XF"
wiesner.data[151:200 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "PX"
wiesner.data[201:250 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "LP"
wiesner.data[251:300 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "PH"

# import SMX measurements
wiesner.data.SMX <- read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/measurements_SMX.csv",
                          skip = 2)
wiesner.data.SMX$replicate <- as.factor(wiesner.data.SMX$replicate)

wiesner.data <- full_join(select(wiesner.data, -technical), select(wiesner.data.SMX, -technical))

wiesner.data$genotype <-
  ordered(
    wiesner.data$genotype,
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
      "cad4xcad5"
    )
  )

# calculate the correct hue on the 360 point circular scale
wiesner.data$hue <- ((wiesner.data$h.stained + 128) / 255 * 360)

wiesner.data$replicate <-
  as.factor(as.character(wiesner.data$replicate))

# calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem
wiesner.data$diff <-
  wiesner.data$OD.stained - wiesner.data$OD.unstained
wiesner.data.bg <- wiesner.data %>%
  filter(cell.type == "PH") %>%
  select(1:3, 9) %>%
  group_by(genotype, replicate) %>%
  summarise(OD.bg = mean(diff, na.rm = TRUE))

wiesner.data.bg$cell.type <- NULL
wiesner.data <-
  merge(
    wiesner.data,
    wiesner.data.bg,
    all = TRUE,
    by = c("genotype", "replicate")
  )
wiesner.data$diff.adj <- wiesner.data$diff - wiesner.data$OD.bg
wiesner.data <- subset(wiesner.data, cell.type != "PH")

# average per replicate (for boxplots)
wiesner.data.pre <-  wiesner.data %>%
  mutate(cell.type = recode(cell.type, "MX" = "PMX"),
         genotype = recode(genotype, "col-0" = "Col-0",
                           "4cl1x2" = "4cl1x4cl2")) %>%
  group_by(genotype, cell.type, replicate) %>%
  summarise(
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  ) %>%
  select(genotype, cell.type, replicate, mean.OD1)

#### import and tidy shape data of the vessels from the Raman measurements ####
raman.irx <- read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/raman_IRX_revisited.csv") %>%
  select(genotype:Perim., Circ., Round, Height, Width) %>%
  filter(object == 2) %>%
  select(-object)
raman.nb <- read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/raman_IRX.csv") %>%
  select(genotype:n_f)
raman.irx <- left_join(raman.irx, raman.nb)
raman.irx$replicate <- as.character(raman.irx$replicate)
raman.irx$technical <- as.character(raman.irx$technical)

#### import and tidy plant height data ####
raman.heights <- read.csv("file:///home/leonard/Documents/Uni/Master/Summer project 16/phenotyping/phenotyping.csv") %>%
  select(genotype, replicate, height) %>%
  mutate(genotype = recode(genotype, "col-0" = "Col-0")) %>%
  rename("Plant.height" = "height") %>%
  rbind(read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/heights_haris.csv"))
raman.heights$replicate <- as.character(raman.heights$replicate)

#### import and tidy Raman data ####
raman.files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/IRX/RAMAN/nuoendagula/",
    pattern = "*.txt",
    recursive = TRUE,
    full.names = TRUE
  )

read_plus <- function(flnm) {
  read_tsv(flnm,
           comment = "#",
           col_names = FALSE,
           skip = 1,
           cols_only(X1 = col_number(),
                     X2 = col_number())
  ) %>%
    mutate(filename = flnm)
}

raman.data <- lapply(raman.files, read_plus) %>%
  bind_rows()

raman.data <- raman.data %>%
  select(c(1:3))

raman.data$filename <- basename(raman.data$filename)

raman.data$filename <-
  str_replace(raman.data$filename, fixed("fah 1"), "fah1")

raman.data <- raman.data %>%
  distinct(X1, filename, .keep_all = TRUE) %>%
  separate(
    filename,
    into = c("genotype", "replicate", "cell.type", "technical"),
    sep = "\\s*#|-|\\s",
    extra = "merge"
  ) %>%
  rename("wavenumber" = X1, "intensity" = X2) %>%
  mutate(
    genotype = recode(
      genotype,
      "4cl1＆2" = "4cl1x4cl2",
      "cad4＆5" = "cad4xcad5",
      "ccoaomt" = "ccoaomt1",
      "ccr1" = "ccr1-3",
      "col.o" = "Col-0",
      "Col.0" = "Col-0",
      "ccr1＆fah1" = "ccr1xfah1"
    ),
    cell.type = recode(cell.type,
                       "px" = "PX",
                       "SX" = "SMX",
                       "MX" = "PMX"),
    replicate = str_extract(replicate,
                            "\\d"),
    technical = str_extract(technical,
                            "(\\d)+")
  ) %>%
  filter(cell.type != "？？")

raman.data$wavenumber <- round(raman.data$wavenumber, digits = 0)

#### baseline correct Raman data ####
raman.data.wide <- raman.data %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate(group_id = row_number()) %>%
  spread(wavenumber, intensity) %>%
  select(-group_id) %>%
  summarise_all(funs(sum(., na.rm = TRUE))) 

raman.data.mat <- as.matrix(raman.data.wide[, -c(1:4)])
rownames(raman.data.mat) <- paste(raman.data.wide$genotype, 
                                           raman.data.wide$cell.type,
                                           raman.data.wide$replicate,
                                           raman.data.wide$technical,
                                           sep = "_")

raman.data.correction <- baseline.als(raman.data.mat)
raman.data.corrected <- data.frame(raman.data.correction$corrected)

raman.data.corrected <- raman.data.corrected %>%
  rownames_to_column() %>%
  separate(rowname, sep = "_", into = c("genotype", "cell.type", "replicate", "technical"), remove = TRUE) %>%
  gather(key = "wavenumber", value = "corrected.intensity", - c(genotype, cell.type, replicate, technical)) %>%
  mutate(wavenumber = str_replace(wavenumber, pattern = fixed("X."), replacement = "-")) %>%
  mutate(wavenumber = str_replace(wavenumber, pattern = fixed("X"), replacement = "")) %>%
  mutate(wavenumber = as.numeric(wavenumber)) %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  inner_join(., raman.data, by = c("genotype", "cell.type", "replicate", "technical", "wavenumber"))

raman.data.corrected$genotype <- ordered(raman.data.corrected$genotype,
                                         levels = c(
                                           "Col-0",
                                           "4cl1",
                                           "4cl2",
                                           "4cl1x4cl2",
                                           "ccoaomt1",
                                           "fah1",
                                           "omt1",
                                           "ccr1-3",
                                           "ccr1xfah1",
                                           "cad4",
                                           "cad5",
                                           "cad4xcad5"
                                         ))

#### average Raman data, calculate integrals and detect peaks ####
raman.data.plot <- raman.data.corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate("1599" = corrected.intensity[wavenumber == 1599],
         "1662" = corrected.intensity[wavenumber == 1662],
         "1119" = corrected.intensity[wavenumber == 1119],
         "1599/1119" = corrected.intensity[wavenumber == 1599] / corrected.intensity[wavenumber == 1119],
         "lig.peak" = MESS::auc(wavenumber, corrected.intensity, from = 1550, to = 1640),
         "cellu.peak" = MESS::auc(wavenumber, corrected.intensity, from = 1110, to = 1130)) %>%
  mutate("1599/1119" = ifelse(`1599/1119` > 30, NA, ifelse(`1599/1119` < -10, NA, `1599/1119`))
  )

raman.data.peaks <- distinct(raman.data.corrected) %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  filter(wavenumber > 1500 & wavenumber < 2000) %>%
  mutate("lig.peak.pos" = wavenumber[which.max(corrected.intensity)]) %>%
  filter(wavenumber == 1599)

raman.data.plot <- inner_join(raman.data.plot, raman.data.peaks) %>%
  ungroup() %>%
  select(-wavenumber, -intensity, -corrected.intensity) %>%
  gather(key = variable, value = value, -c(genotype:technical)) %>%
  group_by(variable) %>%
  mutate(value.scaled = scale_this(value)) %>%
  filter(cell.type != "IF" & genotype != "ccr1xfah1") %>%
  mutate(cell.type = recode(cell.type, "MX" = "PMX"),
         cell.type = ordered(cell.type, levels = c("PX", "PMX", "SMX")))

raman.data.spread <- raman.data.plot %>%
  ungroup() %>%
  group_by(genotype, cell.type, replicate) %>%
  select(-value.scaled) %>%
  spread(variable, value)

raman.data.pre <- raman.data.corrected %>%
  group_by(genotype, cell.type, replicate, wavenumber) %>%
  summarise(
    mean.intensity.pre = mean(corrected.intensity, na.rm = TRUE),
    sd.intensity.pre = sd(corrected.intensity, na.rm = TRUE),
    samples.pre = length(corrected.intensity)
  )

raman.data.avg <- raman.data.pre %>%
  group_by(genotype, cell.type, wavenumber) %>%
  summarise(
    mean.intensity = mean(mean.intensity.pre, na.rm = TRUE),
    sd.intensity = sd(mean.intensity.pre, na.rm = TRUE),
    plants = paste("plants: ", length(mean.intensity.pre)),
    bundles = ifelse(length(mean.intensity.pre) > 1, 
                     paste("cells/plant: ", min(samples.pre), "—", max(samples.pre), sep = ""),
                     paste("cells/plant: ", min(samples.pre)))
  )

#### merge data frames ####
raman.irx <- left_join(raman.irx, raman.data.spread)
raman.irx.switch <- left_join(read_tsv("/home/leonard/Documents/Uni/PhD/IRX/PMX_which_SMX.tsv",
                                       col_types = "cccc"), raman.irx) %>%
  mutate(technical = str_c(technical, "switched", sep = "_"),
         cell.type = "SMX")
raman.irx <- raman.irx %>%
  anti_join(read_tsv("/home/leonard/Documents/Uni/PhD/IRX/PMX_which_SMX.tsv",
                     col_types = "cccc")) %>%
  bind_rows(raman.irx.switch)
raman.irx <- left_join(raman.irx, raman.heights)
raman.irx <- left_join(raman.irx, wiesner.data.pre)
write.csv(raman.irx, file = "SEM_data.csv")

#### create average data frame ####
irx.avg <- read_csv("Wiesner_IRX_data.csv")
irx.avg <- irx.avg %>%
  select(-technical, -number) %>%
  rename(cell.type = object) %>%
  group_by(genotype, replicate, cell.type) %>%
  summarise_all(list(~mean(.)), na.rm = TRUE)
raman.spread.avg <- raman.data.spread %>%
  select(-technical) %>%
  filter(genotype != "ccr1xfah1") %>%
  group_by(genotype, replicate, cell.type) %>%
  summarise_all(list(~mean(.)), na.rm = TRUE)
irx.avg <-left_join(raman.spread.avg, irx.avg)
irx.avg <- left_join(irx.avg, wiesner.data.pre)
write.csv(raman.irx, file = "SEM_data_avg.csv")