library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(showtext)

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
      "ccr1＆fah1" = "ccr1xfah1"
    ),
    cell.type = recode(cell.type,
                       "px" = "PX"),
    replicate = recode(replicate,
                       "＃5" = "5")
  )

raman.data$technical <- str_remove(raman.data$technical, ".txt")

raman.data$genotype <- factor(raman.data$genotype)
raman.data$cell.type <- factor(raman.data$cell.type)
raman.data$replicate <- factor(raman.data$replicate)
raman.data$technical <- factor(raman.data$technical)
raman.data$wavenumber <- round(raman.data$wavenumber, digits = 0)

#### average data ####
raman.data.ratios <- raman.data %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate("rel.cellulose" = intensity / intensity[wavenumber == 381],
         "rel.hemicellulose" = intensity / intensity[wavenumber == 1257],
         "rel.lignin" = intensity / intensity[wavenumber == 1599])

raman.data.pre <- raman.data %>%
  group_by(genotype, cell.type, replicate, wavenumber) %>%
  summarise(
    mean.intensity.pre = mean(intensity, na.rm = TRUE),
    sd.intensity.pre = sd(intensity, na.rm = TRUE)
  )

raman.data.avg <- raman.data.pre %>%
  group_by(genotype, cell.type, wavenumber) %>%
  summarise(
    mean.intensity = mean(mean.intensity.pre, na.rm = TRUE),
    sd.intensity = sd(mean.intensity.pre, na.rm = TRUE),
    samples = paste("n = ", length(mean.intensity.pre))
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
            aes(x = wavenumber, y = mean.intensity),
            size = 0.1) +
  geom_ribbon(
    data = raman.data.pre,
    aes(x = wavenumber,
      fill = replicate,
      ymin = mean.intensity.pre - sd.intensity.pre,
      ymax = mean.intensity.pre + sd.intensity.pre
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
        ymin = mean.intensity - sd.intensity,
        ymax = mean.intensity + sd.intensity
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
    x = 1900,
    y = 19000,
    aes(label = samples),
    stat = "identity",
    family = "Helvetica"
  )

pdf("spectra_RAMAN.pdf", 15, 15)
spectra.WT.MX
dev.off()

rel.lignin <- ggplot(data = subset(raman.data, wavenumber == 1599 & cell.type != "？？"), 
                     aes(x = genotype, y = intensity)) +
  geom_jitter(width = 0.1) +
  geom_violin() +
  # scale_y_continuous(limits = c(-20, 20)) +
  facet_wrap(~ cell.type, ncol = 2)

pdf("lignin_to_cellulose.pdf")
rel.lignin
dev.off()
