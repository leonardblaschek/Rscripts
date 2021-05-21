library(tidyverse)
library(baseline)

#### list Raman spectra files ####
Raman_files <-
  list.files(
    path = "/home/leonard/Downloads/Manuela/",
    pattern = "*.txt",
    recursive = TRUE,
    full.names = TRUE
  )

#### write function to read spectra ####
read_Raman <- function(flnm) {
  read_tsv(flnm,
           comment = "#",
           col_names = c("wavenumber", "intensity"),
           col_types = "cc",
  ) %>%
    mutate(
      filename = basename(flnm)
    ) %>%
    mutate(
      wavenumber = as.numeric(str_replace(wavenumber, fixed(","), fixed("."))),
      intensity = as.numeric(str_replace(intensity, fixed(","), fixed("."))),
    ) %>%
    filter(wavenumber > 300) %>%
    pivot_wider(
      id_cols = filename,
      names_from = wavenumber,
      values_from = intensity
    ) %>%
    group_by(filename) %>%
    nest()
}

Raman_spectra <- map_dfr(Raman_files, read_Raman)

#### baseline correction ####
Raman_als <- function(x) {
  corrected_spectra <- baseline.als(as.matrix(x),
                                    lambda = 5,
                                    p = 0.01,
                                    maxit = 100
  )
  as_tibble(corrected_spectra$corrected, rownames = NA)
}

Raman_spectra_corrected <- Raman_spectra %>%
  mutate(corrected_data = map(data, Raman_als)) %>%
  select(-data) %>%
  unnest() %>%
  pivot_longer(-filename,
               names_to = "wavenumber",
               values_to = "corrected_intensity"
  ) %>%
  drop_na() %>%
  mutate(wavenumber = as.numeric(wavenumber))

#### filtering and scaling ####

# alignment
cellulose_peak <- Raman_spectra_corrected %>%
  filter(round(wavenumber, digits = 0) %in% c(370:390)) %>%
  group_by(filename) %>%
  mutate(
    peak_pos = wavenumber[which.max(corrected_intensity)]
  ) %>%
  select(-wavenumber, -corrected_intensity) %>%
  unique() %>%
  ungroup() %>% 
  mutate(peak_pos_ref = peak_pos[1])

Raman_spectra_aligned <- Raman_spectra_corrected %>%
  left_join(cellulose_peak) %>%
  ungroup() %>%
  mutate(wavenumber = wavenumber + (peak_pos_ref - peak_pos))

# scaling and filtering
Raman_spectra_scaled <- Raman_spectra_aligned %>%
  group_by(filename) %>%
  mutate(
    AUC = MESS::auc(wavenumber, corrected_intensity),
    scaled = corrected_intensity / AUC,
    scaled_lig = corrected_intensity / (
      (corrected_intensity[abs(wavenumber - 1600) == min(abs(wavenumber - 1600))][1] / 2) + corrected_intensity[abs(wavenumber - 1330) == min(abs(wavenumber - 1330))][1]
    )
  )

ggplot(Raman_spectra_scaled,
       aes(x = wavenumber,
           y = scaled,
           group = filename)) +
  geom_line()

Raman_spectra_avg <- Raman_spectra_scaled %>%
  mutate(wavenumber = round(wavenumber, digits = 0)) %>%
  complete(wavenumber = seq(300, 1800)) %>% 
  arrange(filename, wavenumber) %>%
  summarise(
    roll_scaled = slider::slide_index_dbl(scaled, wavenumber, mean, .before = 3, .after = 2),
    roll_scaled_lig = slider::slide_index_dbl(scaled_lig, wavenumber, mean, .before = 2, .after = 2),
    wavenumber = wavenumber
  )