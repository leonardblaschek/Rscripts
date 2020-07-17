library(tidyverse)
library(tukeygrps)

rmn_data_corrected <- read_csv("rmn_data_corrected.csv", col_types = "ccccnnn") %>%
  mutate(
    genotype = ordered(genotype,
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
                       )
    )
  )

rmn_data_corrected <- rmn_data_corrected %>%
  # left_join(lig_peak) %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate(
    AUC = MESS::auc(wavenumber, corrected.intensity),
    scaled = corrected.intensity / AUC,
    scaled_lig = corrected.intensity / (
      (corrected.intensity[wavenumber == 1603] / 2) + corrected.intensity[wavenumber == 1334]
    )
  )

rmn_data_pre <- rmn_data_corrected %>%
  group_by(genotype, cell.type, replicate, wavenumber) %>%
  summarise(
    samples.pre = length(corrected.intensity),
    mean.scaled = mean(scaled, na.rm = TRUE),
    mean.scaled_lig = mean(scaled_lig, na.rm = TRUE),
    # mean.lig_peak = mean(lig_peak, na.rm = TRUE)
  )

rmn_data_avg <- rmn_data_pre %>%
  group_by(genotype, cell.type, wavenumber) %>%
  summarise(
    mean.scaled = mean(mean.scaled, na.rm = TRUE),
    mean.scaled_lig = mean(mean.scaled_lig, na.rm = TRUE),
    # mean.lig_peak = mean(mean.lig_peak, na.rm = TRUE)
  )

#### average Raman data, calculate integrals and detect peaks ####
rmn_sum <- rmn_data_corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  summarise(
    "total_lignin" = (scaled[wavenumber == 1603] / 2) + scaled[wavenumber == 1334],
    "GOH_rel" = scaled_lig[wavenumber == 1664],
    # "GOH.G" = scaled_lig[wavenumber == 1664] / scaled_lig[wavenumber == 1603],
    "GOH_tot" = scaled[wavenumber == 1664],
    "G_tot" = scaled[wavenumber == 1603],
    "G_rel" = scaled_lig[wavenumber == 1603],
    "S.G" = scaled[wavenumber == 1334] / scaled[wavenumber == 1276],
    "GCHO_rel" = scaled_lig[wavenumber == 1625],
    "GCHO_tot" = scaled[wavenumber == 1625],
    "GCHO_1133_rel" = scaled_lig[wavenumber == 1133],
    "GCHO_1133_tot" = scaled[wavenumber == 1133],
    "GCHO_1141_rel" = scaled_lig[wavenumber == 1141],
    "GCHO_1141_tot" = scaled[wavenumber == 1141],
    "S_rel" = scaled_lig[wavenumber == 1334],
    "S_tot" = scaled[wavenumber == 1334],
    # "cellulose" = scaled[wavenumber == 1099] + scaled[wavenumber == 381],
    # "AUC" = mean(AUC)
  ) %>%
  group_by(genotype, cell.type) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

rmn_sd <- rmn_data_corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  summarise(
    "total_lignin" = (scaled[wavenumber == 1603] / 2) + scaled[wavenumber == 1334],
    "GOH_rel" = scaled_lig[wavenumber == 1664],
    # "GOH.G" = scaled_lig[wavenumber == 1664] / scaled_lig[wavenumber == 1603],
    "GOH_tot" = scaled[wavenumber == 1664],
    "G_tot" = scaled[wavenumber == 1603],
    "G_rel" = scaled_lig[wavenumber == 1603],
    "S.G" = scaled[wavenumber == 1334] / scaled[wavenumber == 1276],
    "GCHO_rel" = scaled_lig[wavenumber == 1625],
    "GCHO_tot" = scaled[wavenumber == 1625],
    "GCHO_1133_rel" = scaled_lig[wavenumber == 1133],
    "GCHO_1133_tot" = scaled[wavenumber == 1133],
    "GCHO_1141_rel" = scaled_lig[wavenumber == 1141],
    "GCHO_1141_tot" = scaled[wavenumber == 1141],
    "S_rel" = scaled_lig[wavenumber == 1334],
    "S_tot" = scaled[wavenumber == 1334],
    # "cellulose" = scaled[wavenumber == 1099] + scaled[wavenumber == 381],
    # "AUC" = mean(AUC)
  ) %>%
  group_by(genotype, cell.type) %>%
  summarise_if(is.numeric, sd, na.rm = TRUE)

#### import and tidy Wiesner test data ####
wiesner.data <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
           skip = 2
  )

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
                             skip = 2
)
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
# summary(wiesner.data)

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
wiesner.data.pre <- wiesner.data %>%
  mutate(
    cell.type = recode(cell.type, "MX" = "PMX"),
    genotype = recode(genotype,
                      "col-0" = "Col-0",
                      "4cl1x2" = "4cl1x4cl2"
    )
  ) %>%
  group_by(genotype, cell.type, replicate) %>%
  summarise(
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  ) %>%
  select(genotype, cell.type, replicate, mean.OD1)

wiesner_sum <- wiesner.data.pre %>%
  group_by(genotype, cell.type) %>%
  summarise(wiesner_mean = mean(mean.OD1),
            wiesner_sd = sd(mean.OD1))

phlog_lac <- read_csv("/home/leonard/Documents/Uni/Master/Master thesis/Phenotyping/phlog_lac.csv") %>%
  group_by(genotype, cell.type) %>%
  summarise(wiesner_mean = mean(OD.adj., na.rm = TRUE),
            wiesner_sd = sd(OD.adj., na.rm = TRUE))

wiesner_avg <- wiesner_sum %>%
  bind_rows(phlog_lac %>% 
              filter(genotype == "lac4x17") %>%
              mutate(cell.type = recode(cell.type,
                                        "MX" = "PMX"))) %>%
  select(- wiesner_sd) %>%
  filter(cell.type %in% c("PMX", "PX", "IF", "SMX"))

wiesner_sd <- wiesner_sum %>%
  bind_rows(phlog_lac %>% 
              filter(genotype == "lac4x17") %>%
              mutate(cell.type = recode(cell.type,
                                        "MX" = "PMX"))) %>%
  select(- wiesner_mean) %>%
  filter(cell.type %in% c("PMX", "PX", "IF", "SMX"))

thio_mean <- full_join(rmn_sum, wiesner_avg) 

thio_sd <- full_join(rmn_sd, wiesner_sd) 

write_csv(thio_mean, path = "histo_mean.csv")

write_csv(thio_sd, path = "histo_sd.csv")

rmn_pre <- rmn_data_corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  summarise(
    "GCHO_tot" = scaled[wavenumber == 1625]
  ) %>%
  group_by(genotype, cell.type, replicate) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

wiesner_ratio <- wiesner.data.pre %>%
  left_join(rmn_pre) %>%
  filter(cell.type %in% c("PMX", "PX", "SMX", "IF")) %>%
  group_by(genotype, cell.type, replicate) %>%
  mutate(ratio = GCHO_tot / mean.OD1) 
# %>%
#   group_by(genotype, cell.type) %>%
#   summarise(mean_ratio = mean(ratio, na.rm = T),
#             sd_ratio = sd(ratio, na.rm = T))

write_csv(wiesner_ratio, path = "wiesner_ratio.csv")

wiesner_stats <- wiesner_ratio %>%
  unite("geno_cell", genotype, cell.type) %>%
  letter_groups(.,
                               ratio,
                               geno_cell,
                               "tukey")


#### spectra ####

planta <- rmn_data_avg %>%
  filter(genotype %in% c("Col-0", "cad4xcad5"))
vitro <- rmn_models %>% filter(gamma == "cho")

write_csv(planta, path = "planta.csv")
write_csv(vitro, path = "vitro.csv")