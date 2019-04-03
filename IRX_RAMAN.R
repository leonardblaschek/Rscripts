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

###############################
# functions for scaling and statistics
###############################
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

#### import and tidy data ####
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
                       "SX" = "SMX"),
    replicate = str_extract(replicate,
                       "\\d"),
    technical = str_extract(technical,
                            "(\\d)+")
  ) %>%
  filter(cell.type != "？？")



# raman.data$genotype <- factor(raman.data$genotype)
# raman.data$cell.type <- factor(raman.data$cell.type)
# raman.data$replicate <- factor(raman.data$replicate)
# raman.data$technical <- factor(raman.data$technical)
raman.data$wavenumber <- round(raman.data$wavenumber, digits = 0)

#### baseline correct ####
raman.data.corrected <- raman.data %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate(group_id = row_number()) %>%
  spread(wavenumber, intensity) %>%
  select(-group_id) %>%
  summarise_all(funs(sum(., na.rm = TRUE))) 

raman.data.corrected.mx <- as.matrix(raman.data.corrected[, -c(1:4)])
rownames(raman.data.corrected.mx) <- paste(raman.data.corrected$genotype, 
                                        raman.data.corrected$cell.type,
                                        raman.data.corrected$replicate,
                                        raman.data.corrected$technical,
                                        sep = "_")

corrected.spectra <- baseline.als(raman.data.corrected.mx)
corrected <- data.frame(corrected.spectra$corrected)

raman.data.corrected <- corrected %>%
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

#### average data ####
raman.data.ratios <- raman.data.corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate("rel.381" = corrected.intensity / corrected.intensity[wavenumber == 381],
         "rel.1119" = corrected.intensity / corrected.intensity[wavenumber == 1119],
         "rel.1599" = corrected.intensity / corrected.intensity[wavenumber == 1599],
         "rel.1257" = corrected.intensity / corrected.intensity[wavenumber == 1257],
         "lig.peak" = MESS::auc(wavenumber, corrected.intensity, from = 1550, to = 1640))

raman.data.plot <- raman.data.corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate("1599 ~ cm ^ -1" = corrected.intensity[wavenumber == 1599],
         "1662 ~ cm ^ -1" = corrected.intensity[wavenumber == 1662],
         "1119 ~ cm ^ -1" = corrected.intensity[wavenumber == 1119],
         "1599 ~ cm ^ -1/1119 ~cm ^ -1" = corrected.intensity[wavenumber == 1599] / corrected.intensity[wavenumber == 1119],
         "Integrated~lignin~peak" = MESS::auc(wavenumber, corrected.intensity, from = 1550, to = 1640),
         "Integrated~cellulose~peak" = MESS::auc(wavenumber, corrected.intensity, from = 1110, to = 1130)) %>%
  mutate("1599 ~ cm ^ -1/1119 ~cm ^ -1" = ifelse(`1599 ~ cm ^ -1/1119 ~cm ^ -1` > 30, NA, ifelse(`1599 ~ cm ^ -1/1119 ~cm ^ -1` < -10, NA, `1599 ~ cm ^ -1/1119 ~cm ^ -1`))
             )

raman.data.peaks <- raman.data.corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  filter(wavenumber > 1500 & wavenumber < 2000) %>%
  mutate("Lignin~peak~position" = wavenumber[which.max(corrected.intensity)]) %>%
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

raman.letters <- raman.data.plot %>%
  group_by(cell.type, variable) %>%
  do(data.frame(tukey(.)))

raman.letters$value <- ifelse(raman.letters$variable == "Lignin~peak~position", 1585, raman.letters$value)


raman.summary <- ggplot(data = raman.data.plot, aes(x = genotype, y = value)) +
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
  theme(strip.text = element_text(face = "italic")) +
  facet_grid(variable ~ cell.type,
             scales = "free_y",
             labeller = label_parsed)
  
pdf("raman_summary.pdf", width = 10, height = 12)
raman.summary
dev.off()
