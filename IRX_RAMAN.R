library(baseline)
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


#### import Helvetica Neue ####
font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

#### import and tidy data ####
raman.files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/IRX/RAMAN",
    pattern = "*.txt",
    recursive = TRUE,
    full.names = TRUE
  )

read_plus <- function(flnm) {
  read_tsv(flnm,
           comment = "#",
           col_names = FALSE,
           skip = 1) %>%
    mutate(filename = flnm)
}

raman.data <- lapply(raman.files, read_plus) %>%
  bind_rows()

raman.data <- raman.data %>%
  select(c(1:3)) %>%
  filter(!str_detect(filename, "baseline corrected"))

raman.data$filename <- basename(raman.data$filename)

raman.data$filename <-
  str_replace(raman.data$filename, fixed("fah 1"), "fah1")

raman.data <- raman.data %>%
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

#### average data ####
raman.data.ratios <- raman.data.corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate("rel.381" = corrected.intensity / corrected.intensity[wavenumber == 381],
         "rel.1119" = corrected.intensity / corrected.intensity[wavenumber == 1119],
         "rel.1599" = corrected.intensity / corrected.intensity[wavenumber == 1599])

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
                     paste("bundles: ", min(samples.pre), "—", max(samples.pre), sep = ""),
                     paste("bundles: ", min(samples.pre)))
  )

spectra.WT.MX <- ggplot() +
  geom_vline(xintercept = 382, size = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 1600, size = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 1256, size = 0.2, alpha = 0.5) +
  annotate("text", x = 382, y = 5000, 
           label = "Cellulose~(382~cm^-1)", 
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
  theme_few() +
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

rel.lignin <- ggplot(data = subset(raman.data.corrected, wavenumber == 1119),
                     aes(x = genotype, y = log(corrected.intensity), fill = replicate)) +
  geom_jitter(width = 0.25, shape = 21, stroke = 0.1, size = 2) +
  # geom_violin(draw_quantiles = 0.5, fill = NA) +
  geom_boxplot(fill = NA) +
  # scale_y_continuous(limits = c(-50, 50)) +
  scale_fill_viridis_d() +
  facet_wrap(~ cell.type, ncol = 2) +
  theme_few() +
  theme(text = element_text(family = "Helvetica"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf("lignin_to_cellulose.pdf")
rel.lignin
dev.off()
