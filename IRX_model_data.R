library(baseline)
library(agricolae)
library(ggthemes)
library(showtext)
library(zoo)
library(cowplot)
library(tidyverse)


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

wiesner.data <- rbind(select(wiesner.data, -technical), select(wiesner.data.SMX, -technical))

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
summary(wiesner.data)

wiesner.data.bg <- wiesner.data %>%
  filter(cell.type == "PH") %>%
  select(1:2, 9) %>%
  ungroup() %>%
  group_by(genotype, replicate) %>%
  summarise(OD.bg = mean(diff, na.rm = TRUE))

wiesner.data <-
  full_join(
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
    SD.OD1 = sd(diff.adj, na.rm = TRUE)) %>%
  mutate(value.scaled = scale_this(mean.OD1)) %>%
  select(genotype, cell.type, replicate, mean.OD1, value.scaled)

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
write_csv(raman.irx, "raman_irx_shape.csv")

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
write_csv(raman.data.corrected, "raman_spectra_corrected.csv")

#### average Raman data, calculate integrals and detect peaks ####
raman.data.plot <- raman.data.corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate("1603" = corrected.intensity[wavenumber == 1603],
         "1662" = corrected.intensity[wavenumber == 1662],
         "1119" = corrected.intensity[wavenumber == 1119],
         "1603/1119" = corrected.intensity[wavenumber == 1603] / corrected.intensity[wavenumber == 1119],
         "1621/1603" = corrected.intensity[wavenumber == 1621] / corrected.intensity[wavenumber == 1603],
         "1340/1603" = corrected.intensity[wavenumber == 1340] / corrected.intensity[wavenumber == 1603],
         "1670/1603" = corrected.intensity[wavenumber == 1670] / corrected.intensity[wavenumber == 1603],
         "lig.peak" = MESS::auc(wavenumber, corrected.intensity, from = 1550, to = 1640),
         "cellu.peak" = MESS::auc(wavenumber, corrected.intensity, from = 1110, to = 1130))

raman.data.peaks <- distinct(raman.data.corrected) %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  filter(wavenumber > 1500 & wavenumber < 2000) %>%
  mutate("lig.peak.pos" = wavenumber[which.max(corrected.intensity)]) %>%
  filter(wavenumber == 1603)

raman.data.plot <- inner_join(raman.data.plot, raman.data.peaks) %>%
  ungroup() %>%
  select(-wavenumber, -intensity, -corrected.intensity) %>%
  gather(key = variable, value = value, -c(genotype:technical)) %>%
  group_by(variable) %>%
  mutate(value.scaled = scale_this(value)) %>%
  filter(
    # cell.type != "IF" & 
      genotype != "ccr1xfah1"
    ) %>%
  mutate(cell.type = recode(cell.type, "MX" = "PMX"),
         cell.type = ordered(cell.type, levels = c("IF", "PX", "PMX", "SMX")))

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

#### plotting ####

# spectra overview
spectra.WT.MX <- ggplot() +
  geom_vline(xintercept = 1120, size = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 1600, size = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 1256, size = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 1550, size = 0.1, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 1640, size = 0.1, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 1110, size = 0.1, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 1130, size = 0.1, alpha = 0.5, linetype = 2) +
  annotate("text", x = 1120, y = 5000, 
           label = "Cellulose~(1120~cm^-1)", 
           parse = TRUE,
           family = "Helvetica",
           size = 2,
           hjust = 0) +
  annotate("text", x = 1600, y = 15000, 
           label = "Lignin~(1600~cm^-1)", 
           parse = TRUE,
           family = "Helvetica",
           size = 2,
           hjust = 0) +
  annotate("text", x = 1256, y = -1000, 
           label = "Hemicellulose~(1256~cm^-1)", 
           parse = TRUE,
           family = "Helvetica",
           size = 2,
           hjust = 0) +
  geom_line(data = raman.data.avg,
            aes(x = wavenumber, y = rollmean(mean.intensity, 5, na.pad=TRUE)),
            size = 0.1) +
  geom_ribbon(
    data = raman.data.pre,
    aes(x = wavenumber,
        fill = replicate,
        ymin = rollmean(mean.intensity.pre, 5, na.pad=TRUE) - rollmean(sd.intensity.pre, 5, na.pad=TRUE),
        ymax = rollmean(mean.intensity.pre, 5, na.pad=TRUE) + rollmean(sd.intensity.pre, 5, na.pad=TRUE)
    ),
    alpha = 0.25
  ) +
  geom_ribbon(
    data = raman.data.avg,
    fill = NA,
    colour = "black",
    linetype = 2,
    size = 0.05,
    aes(x = wavenumber,
        ymin = rollmean(mean.intensity, 5, na.pad=TRUE) - rollmean(sd.intensity, 5, na.pad=TRUE),
        ymax = rollmean(mean.intensity, 5, na.pad=TRUE) + rollmean(sd.intensity, 5, na.pad=TRUE)
    ),
    alpha = 0.25
  ) +
  facet_grid(genotype ~ cell.type) +
  scale_x_continuous(limits = c(300, 2000)) +
  scale_y_continuous(limits = c(-1000, 20000)) +
  scale_fill_viridis_d() +
  theme_leo() +
  theme(text = element_text(family = "Helvetica")) +
  geom_text(
    data = subset(raman.data.avg, wavenumber == 1430),
    x = 2050,
    y = 19500,
    aes(label = plants),
    stat = "identity",
    family = "Helvetica",
    size = 2.5,
    hjust = 1
  ) +
  geom_text(
    data = subset(raman.data.avg, wavenumber == 1430),
    x = 2050,
    y = 17500,
    aes(label = bundles),
    stat = "identity",
    family = "Helvetica",
    size = 2.5,
    hjust = 1
  )

pdf("spectra_RAMAN.pdf", 15, 15)
spectra.WT.MX
dev.off()

# raman spectra genoype comparison
spectra.comp <- ggplot(data = filter(raman.data.avg, genotype == "Col-0" & cell.type == "PMX")) +
  geom_vline(xintercept = 1120, size = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 1600, size = 0.2, alpha = 0.5) +
  # geom_vline(xintercept = 1256, size = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 1550, size = 0.1, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 1640, size = 0.1, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 1110, size = 0.1, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 1130, size = 0.1, alpha = 0.5, linetype = 2) +
  annotate("text", x = 1120, y = 5000,
           label = "Cellulose~(1120~cm^-1)",
           parse = TRUE,
           family = "Helvetica",
           size = 4,
           hjust = 0) +
  annotate("text", x = 1600, y = 11000,
           label = "Lignin~(1600~cm^-1)",
           parse = TRUE,
           family = "Helvetica",
           size = 4,
           hjust = 0) +
  # annotate("text", x = 1256, y = -1000, 
  #          label = "Hemicellulose~(1256~cm^-1)", 
  #          parse = TRUE,
  #          family = "Helvetica",
  #          size = 2,
  #          hjust = 0) +
  geom_line(aes(x = wavenumber,
                y = rollmean(mean.intensity, 5, na.pad=TRUE)),
            size = 0.3) +
  # geom_ribbon(
  #   data = raman.data.pre,
  #   aes(x = wavenumber,
  #       fill = replicate,
  #       ymin = rollmean(mean.intensity.pre, 5, na.pad=TRUE) - rollmean(sd.intensity.pre, 5, na.pad=TRUE),
  #       ymax = rollmean(mean.intensity.pre, 5, na.pad=TRUE) + rollmean(sd.intensity.pre, 5, na.pad=TRUE)
  #   ),
  #   alpha = 0.25
  # ) +
  geom_ribbon(
    # colour = "black",
    # linetype = 2,
    # size = 0.05,
    aes(x = wavenumber,
        ymin = rollmean(mean.intensity, 5, na.pad=TRUE) - rollmean(sd.intensity, 5, na.pad=TRUE),
        ymax = rollmean(mean.intensity, 5, na.pad=TRUE) + rollmean(sd.intensity, 5, na.pad=TRUE)
    ),
    alpha = 0.25
  ) +
  facet_wrap(~ cell.type, nrow = 1) +
  scale_x_continuous(limits = c(300, 2000)) +
  scale_y_continuous(limits = c(-1000, 11000)) +
  # scale_fill_viridis_d() +
  # scale_fill_brewer(palette = "Set1") +
  # scale_colour_brewer(palette = "Set1") +
  theme_leo() +
  labs(x = bquote(Wavenumber~(cm^-1))) +
  theme(text = element_text(family = "Helvetica"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = c(0.05, 0.85),
        legend.key.width = unit(3, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.title = element_blank(),
        axis.title = element_text(),
        axis.title.y = element_blank())
  # geom_text(
  #   data = subset(raman.data.avg, wavenumber == 1430),
  #   x = 2050,
  #   y = 19500,
  #   aes(label = plants),
  #   stat = "identity",
  #   family = "Helvetica",
  #   size = 2.5,
  #   hjust = 1
  # ) +
  # geom_text(
  #   data = subset(raman.data.avg, wavenumber == 1430),
  #   x = 2050,
  #   y = 17500,
  #   aes(label = bundles),
  #   stat = "identity",
  #   family = "Helvetica",
  #   size = 2.5,
  #   hjust = 1
  # )

pdf("comp_spectra_RAMAN.pdf", 10, 3)
spectra.comp
dev.off()

# raman data overview
raman.letters <- filter(raman.data.plot, variable %in% c("1603", "1119", "1603/1119", "1340/1603")) %>%
  group_by(cell.type, variable) %>%
  do(data.frame(tukey(.)))

raman.letters$value <- ifelse(raman.letters$variable == "lig.peak.pos", 1585, raman.letters$value)


raman.summary <- ggplot(data = filter(raman.data.plot, variable %in% c("1603", "1119", "1603/1119", "1340/1603")), aes(x = genotype, y = value)) +
  geom_jitter(
    aes(fill = value.scaled),
    shape = 21,
    width = 0.1,
    alpha = 0.9,
    size = 2,
    stroke = 0.25
  ) +
  # geom_violin(draw_quantiles = 0.5, adjust = 1.5, fill = rgb(1,1,1,0.5)) +
  geom_boxplot(fill = rgb(1,1,1,0.5), outlier.alpha = 0) +
  geom_text(data = raman.letters,
            aes(label = groups),
            angle = 90,
            hjust = 1,
            family = "Helvetica") +
  scale_fill_distiller(palette = "RdBu", name = "Z-score by\nrow", limits = c(-10, 10)) +
  scale_y_continuous(expand = expand_scale(mult = c(0.2,0.05))) +
  scale_x_discrete(
    labels = c(
      "Col-0",
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(paste(italic("4cl1"), "x", italic("4cl2"))),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(italic("cad4")),
      expression(italic("cad5")),
      expression(paste(italic("cad4"), "x", italic("cad5")))
    )
  ) +
  theme_leo() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none") +
  facet_grid(variable ~ cell.type,
             scales = "free_y")

pdf("raman_summary.pdf", width = 10, height = 6)
raman.summary
dev.off()

# wiesner data overview
# 
# wiesner.letters <- filter(wiesner.data.pre, cell.type %in% c("PX", "PMX", "SMX")) %>%
#   group_by(cell.type) %>%
#   do(data.frame(tukey(.)))
wiesner.data.pre$cell.type <- ordered(wiesner.data.pre$cell.type, levels = c("IF", "LP", "XF", "PX", "PMX", "SMX"))
wiesner.summary <- ggplot(data = filter(wiesner.data.pre, cell.type %in% c("PX", "PMX", "SMX") & genotype != "ccr1xfah1"), aes(x = genotype, y = mean.OD1)) +
  geom_jitter(
    aes(fill = value.scaled),
    shape = 21,
    width = 0.1,
    alpha = 0.9,
    size = 2,
    stroke = 0.25
  ) +
  # geom_violin(draw_quantiles = 0.5, adjust = 1.5, fill = rgb(1,1,1,0.5)) +
  geom_boxplot(fill = rgb(1,1,1,0.5), outlier.alpha = 0) +
  # geom_text(data = wiesner.letters,
  #           aes(label = groups),
  #           angle = 90,
  #           hjust = 1,
  #           family = "Helvetica") +
  scale_fill_distiller(palette = "RdBu", name = "Z-score by\nrow", limits = c(-6, 6)) +
  scale_y_continuous(expand = expand_scale(mult = c(0.2,0.05))) +
  scale_x_discrete(
    labels = c(
      "Col-0",
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(paste(italic("4cl1"), "x", italic("4cl2"))),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(italic("cad4")),
      expression(italic("cad5")),
      expression(paste(italic("cad4"), "x", italic("cad5")))
    )
  ) +
  theme_leo() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none") +
  facet_grid(~ cell.type)

pdf("wiesner_summary.pdf", width = 10, height = 3)
wiesner.summary
dev.off()

raman.irx$genotype <- ordered(raman.irx$genotype,
                              levels = c(
                                "Col-0",
                                "4cl1",
                                "4cl2",
                                "4cl1x4cl2",
                                "ccoaomt1",
                                "fah1",
                                "omt1",
                                "ccr1-3",
                                "cad4",
                                "cad5",
                                "cad4xcad5"
                              ))

raman.irx$cell.type <- ordered(raman.irx$cell.type,
                               levels = c(
                                 "PX",
                                 "PMX",
                                 "SMX"
                               ))

# irx neighbours
irx.neighbours <-
  ggplot(data = select(raman.irx, genotype:cell.type, n_f, n_v, n_p) %>% 
           gather(key = "adjacency", value = "proportion",
                   n_f, n_v, n_p) %>%
           group_by(genotype, cell.type, adjacency) %>%
           summarise(proportion.mean = mean(proportion, na.rm = TRUE),
                     proportion.n = n()),
         aes(x = genotype, y = proportion.mean)) +
  geom_bar(
    aes(fill = adjacency),
    stat = "identity",
    # width = 0.1,
    alpha = 0.9,
    size = 2
    # stroke = 0.25
  ) +
  scale_fill_few() +
  # scale_fill_distiller(palette = "RdBu", name = "Z-score by\nrow") +
  # scale_y_continuous(expand = expand_scale(mult = c(0.28,0.05))) +
  scale_x_discrete(
    labels = c(
      "Col-0",
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(paste(italic("4cl1"), "x", italic("4cl2"))),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(italic("cad4")),
      expression(italic("cad5")),
      expression(paste(italic("cad4"), "x", italic("cad5")))
    )
  ) +
  geom_text(aes(x = genotype, y = 0, label = proportion.n)) +
  theme_leo() +
  theme(legend.position = "bottom") +
  facet_wrap(~ cell.type, nrow = 1) 
# coord_flip()

pdf("irx_overview_5.pdf", width = 10, height = 6)
irx.neighbours
dev.off()